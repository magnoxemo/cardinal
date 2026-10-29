#ifdef ENABLE_OPENMC_COUPLING
#include "FETMeshTally.h"

#include "AuxiliarySystem.h"

#include "openmc/fet/dubiner.h"

registerMooseObject("CardinalApp", FETMeshTally);

InputParameters
FETMeshTally::validParams()
{
  auto params = MeshTally::validParams();
  params.addClassDescription(
      "A class which implements unstructured mesh tallies expanded in a functional expansion "
      "basis within each element. Only linear tetrahedral elements are supported.");
  params.addRequiredParam<int>("order", "Functional expansion order used on every element");
  return params;
}

FETMeshTally::FETMeshTally(const InputParameters & parameters)
  : MeshTally(parameters), _order(getParam<int>("order"))
{
  if (_order < 0 || _order > openmc::fet::MAX_ORDER)
    paramError("order",
               "The expansion order must be between 0 and " +
                   std::to_string(openmc::fet::MAX_ORDER) + "!");

  // The filter has one bin per (element, mode), so the element-wise copy and relaxation
  // machinery of the mesh tally does not apply.
  if (_use_dof_map)
    mooseError("Functional expansion mesh tallies do not yet support adaptivity or block "
               "restriction!");
  if (_instance != 0)
    mooseError("Functional expansion mesh tallies do not support multiple mesh translations!");
}

openmc::MeshFilter *
FETMeshTally::createMeshFilter()
{
  _mesh_fet_filter = dynamic_cast<openmc::MeshFETFilter *>(openmc::Filter::create("meshfet"));
  return _mesh_fet_filter;
}

void
FETMeshTally::configureMeshFilter()
{
  // must happen before the tally allocates its result
  // I am leaving that as a future note to myself that doing uniform P refinement
  // is one of the dumbest thing I can do. But that's a worry for post prelim Walid
  _mesh_fet_filter->set_order(_order);
}

Real
FETMeshTally::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                unsigned int local_score,
                                const std::vector<OMCTensor> & tally_vals,
                                bool norm_by_src_rate)
{
  Real total = 0.0;

  const auto n_bins = _mesh_fet_filter->n_bins();
  const auto n_elem = _mesh_fet_filter->n_elements();

  Real factor = 1.0;
  if (norm_by_src_rate)
  {
    factor = _openmc_problem.tallyMultiplier(_tally_score[local_score],
                                             _local_mean_tally[local_score]) *
             _openmc_problem.scaling() * _openmc_problem.scaling() * _openmc_problem.scaling();

    _modal_values.resize(_tally_score.size());
    _modal_values[local_score].assign(_num_ext_filter_bins, std::vector<Real>(n_bins, 0.0));
  }

  for (unsigned int ext_bin = 0; ext_bin < _num_ext_filter_bins; ++ext_bin)
  {
    for (int e = 0; e < n_elem; ++e)
    {
      const auto first = _mesh_fet_filter->offset(e); // the first index of the element offset
      const auto n_modes = openmc::fet::n_modes_tet(_mesh_fet_filter->order(e));
      if (n_modes == 0)
        continue;

      std::vector<Real> values(n_modes);
      for (int i = 0; i < n_modes; ++i)
        values[i] = tally_vals[local_score](ext_bin * n_bins + first + i) * factor;

      total += _ext_bins_to_skip[ext_bin] ? 0.0 : tally_vals[local_score](ext_bin * n_bins + first);

      auto var = var_numbers[_num_ext_filter_bins * local_score + ext_bin];

      // the variable is a discontinuous linear Lagrange variable, so there is one dof per
      // vertex; evaluate the expansion at each vertex (the filter removes its translation
      // internally, so add it to the position, same as in reconstruct())
      const auto * elem = _openmc_problem.getMooseMesh().queryElemPtr(e);
      if (_openmc_problem.isLocalElem(elem))
      {
        auto & solution = _aux.solution();
        for (unsigned int n = 0; n < elem->n_vertices(); ++n)
        {
          const auto r = elem->point(n) * _openmc_problem.scaling() + _mesh_translation;
          const Real nodal = _mesh_fet_filter->reconstruct(
              e, openmc::Position(r(0), r(1), r(2)), values.data());
          solution.set(elem->dof_number(_aux.number(), var, n), nodal);
        }
      }

      if (norm_by_src_rate)
        std::copy(values.begin(), values.end(), _modal_values[local_score][ext_bin].begin() + first);
    }
  }

  return total;
}

Real
FETMeshTally::reconstruct(const libMesh::Elem * elem,
                          const libMesh::Point & p,
                          unsigned int local_score,
                          unsigned int ext_bin) const
{
  const int e = elem->id();
  if (local_score >= _modal_values.size() || ext_bin >= _modal_values[local_score].size())
    mooseError("The tally results have not been stored yet; cannot reconstruct the expansion.");
  if (_mesh_fet_filter->order(e) < 0)
    return 0.0;

  // positions in OpenMC are scaled; the filter removes its own translation internally, so add
  // it here (the same convention as checkMeshTemplateAndTranslations)
  const auto r = p * _openmc_problem.scaling() + _mesh_translation;
  return _mesh_fet_filter->reconstruct(
      e,
      openmc::Position(r(0), r(1), r(2)),
      _modal_values[local_score][ext_bin].data() + _mesh_fet_filter->offset(e));
}

#endif
