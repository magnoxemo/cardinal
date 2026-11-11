#ifdef ENABLE_OPENMC_COUPLING
#include "FETSurfaceTally.h"

#include <array>
#include <cmath>
#include <map>
#include <set>

#include "AuxiliarySystem.h"
#include "MooseMesh.h"
#include "libmesh/mesh_tools.h"

#include "openmc/fet/fet_element.h"
#include "openmc/surface.h"


#ifdef OPENMC_XDG_ENABLED
#include "openmc/xdg.h"
#include "xdg/xdg.h"
#endif

registerMooseObject("CardinalApp", FETSurfaceTally);

InputParameters
FETSurfaceTally::validParams()
{
  auto params = TallyBase::validParams();
  params.addClassDescription(
      "A class which implements surface tallies on XDG surface-mesh surfaces, expanded in a "
      "functional expansion basis on each (triangular) facet.");
  params.addRequiredParam<std::vector<int32_t>>(
      "surface_ids", "OpenMC IDs of the XDG surface-mesh surfaces to tally on");
  params.addRequiredParam<int>("order", "Functional expansion order used on every facet");
  return params;
}

FETSurfaceTally::FETSurfaceTally(const InputParameters & parameters)
  : TallyBase(parameters),
    _surface_ids(getParam<std::vector<int32_t>>("surface_ids")),
    _order(getParam<int>("order"))
{
  if (_surface_ids.empty())
    paramError("surface_ids", "At least one surface must be provided!");

  if (_order < 0 || _order > openmc::fet::MAX_ORDER)
    paramError("order",
               "The expansion order must be between 0 and " +
                   std::to_string(openmc::fet::MAX_ORDER) + "!");

  // A surface tally is not a volumetric source, so there is no unfiltered global tally to
  // normalize against or to check the sum with.
  if (_needs_global_tally)
    mooseError("Surface FET tallies cannot be normalized by, or checked against, a global tally; "
               "set 'normalize_by_global_tally = false' and 'check_tally_sum = false'!");

  if (_is_adaptive)
    mooseError("Surface FET tallies do not support adaptivity!");

  // OpenMC only scores 'flux' (surface-crossing flux) or 'current' on a surface filter
  for (const auto & score : _tally_score)
    if (score != "flux")
      paramError("score", "Surface FET tallies only support the 'flux' score!");

  // Surface crossings are scored by the analog estimator.
  if (isParamValid("estimator"))
  {
    if (_estimator != openmc::TallyEstimator::ANALOG)
      paramError("estimator", "Surface FET tallies require the analog estimator!");
  }
  else
    _estimator = openmc::TallyEstimator::ANALOG;
}

void
FETSurfaceTally::initializeTally()
{
  TallyBase::initializeTally();

  // The scores were set before the filters, so OpenMC did not see the surface filter and left
  // the tally as a volume tally; set them again now that the filters are in place.
  _local_tally->set_scores(_tally_score);
  if (_local_tally->type_ != openmc::TallyType::SURFACE)
    mooseError("Internal error: the surface FET tally was not classified as a surface tally by "
               "OpenMC.");
}

std::pair<unsigned int, openmc::Filter *>
FETSurfaceTally::spatialFilter()
{
  _surface_fet_filter =
      dynamic_cast<openmc::SurfaceFETFilter *>(openmc::Filter::create("surfacefet"));
  if (!_surface_fet_filter)
    mooseError("Internal error: could not create a surface FET filter.");

  // Convert the OpenMC surface IDs to indices into model::surfaces
  std::vector<int32_t> indices;
  for (const auto id : _surface_ids)
  {
    auto search = openmc::model::surface_map.find(id);
    if (search == openmc::model::surface_map.end())
      paramError("surface_ids", "Surface " + std::to_string(id) + " does not exist in OpenMC!");
    indices.push_back(search->second);
  }

  // The order must be set before the facets are built, so that it applies to all of them.
  _surface_fet_filter->set_order(_order);
  _surface_fet_filter->set_surfaces(openmc::span<int32_t>(indices.data(), indices.size()));

  buildFacetSideMap(indices);

  return std::make_pair(openmc::model::tally_filters.size() - 1, _surface_fet_filter);
}

void
FETSurfaceTally::buildFacetSideMap(const std::vector<int32_t> & surface_indices)
{
#ifdef OPENMC_XDG_ENABLED
  const auto & moose_mesh = _openmc_problem.getMooseMesh();
  const Real scaling = _openmc_problem.scaling();

  // XDG reads the mesh through its own libMesh::Mesh, which may renumber elements and nodes, so
  // the facets are matched to the sides of the [Mesh] by geometry: the centroid of a facet (in
  // cm) is the centroid of a side of the [Mesh] times the scaling. Bin the side centroids in a
  // grid and look in the neighboring bins for a match.
  libMesh::BoundingBox bbox = libMesh::MeshTools::create_bounding_box(moose_mesh.getMesh());
  const Real extent = std::max((bbox.max() - bbox.min()).norm() * scaling, 1e-12);
  const Real cell = 1e-3 * extent;
  const Real tol = 1e-6 * extent;

  using Key = std::array<long long, 3>;
  auto key_of = [cell](const Point & p)
  {
    return Key{static_cast<long long>(std::floor(p(0) / cell)),
               static_cast<long long>(std::floor(p(1) / cell)),
               static_cast<long long>(std::floor(p(2) / cell))};
  };

  struct SideCentroid
  {
    Point centroid;
    FacetSide side;
  };
  std::map<Key, std::vector<SideCentroid>> grid;
  for (const auto * elem : moose_mesh.getMesh().active_element_ptr_range())
    for (unsigned int s = 0; s < elem->n_sides(); ++s)
    {
      const auto centroid = elem->side_ptr(s)->vertex_average() * scaling;
      grid[key_of(centroid)].push_back(
          {centroid, {elem->id(), static_cast<unsigned short>(s)}});
    }

  _facet_sides.clear();
  _facet_sides.reserve(_surface_fet_filter->n_elements());

  // The facets of the filter are numbered surface after surface, in the order of
  // get_surface_faces() of the XDG mesh manager; walk them in the same order.
  for (const auto index : surface_indices)
  {
    const auto * surface =
        dynamic_cast<const openmc::XDGSurface *>(openmc::model::surfaces[index].get());
    if (!surface)
      mooseError("Surface ", openmc::model::surfaces[index]->id_, " is not an XDG surface!");

    const auto & mm = *surface->xdg_ptr()->mesh_manager();
    for (const auto face : mm.get_surface_faces(surface->xdg_id()))
    {
      const auto vertices = mm.face_vertex_coordinates(face);
      Point facet_centroid;
      for (const auto & v : vertices)
        facet_centroid += Point(v.x, v.y, v.z);
      facet_centroid /= static_cast<Real>(vertices.size());

      const auto k = key_of(facet_centroid);
      const SideCentroid * best = nullptr;
      Real best_distance = tol;
      for (long long i = -1; i <= 1; ++i)
        for (long long j = -1; j <= 1; ++j)
          for (long long l = -1; l <= 1; ++l)
          {
            auto it = grid.find(Key{k[0] + i, k[1] + j, k[2] + l});
            if (it == grid.end())
              continue;
            for (const auto & candidate : it->second)
            {
              const Real d = (candidate.centroid - facet_centroid).norm();
              if (d <= best_distance)
              {
                best = &candidate;
                best_distance = d;
              }
            }
          }

      if (!best)
        mooseError("Facet ",
                   face,
                   " of surface ",
                   surface->id_,
                   " (centroid ",
                   facet_centroid,
                   " cm) does not coincide with a side of the [Mesh]; the [Mesh] must have the "
                   "same geometry as the mesh OpenMC loaded!");

      _facet_sides.push_back(best->side);
    }
  }

  if (static_cast<int>(_facet_sides.size()) != _surface_fet_filter->n_elements())
    mooseError("Internal error: found ",
               _facet_sides.size(),
               " facet sides for ",
               _surface_fet_filter->n_elements(),
               " facets.");
#else
  mooseError("Surface FET tallies require OpenMC to be built with XDG (OPENMC_USE_XDG=ON)!");
#endif
}
Real
FETSurfaceTally::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                   unsigned int local_score,
                                   const std::vector<OMCTensor> & tally_vals,
                                   bool norm_by_src_rate)
{
  Real total = 0.0;

  const auto n_bins = _surface_fet_filter->n_bins();
  const auto n_elem = _surface_fet_filter->n_elements();

  // per unit area; lengths in OpenMC are scaled
  Real factor = 1.0;
  if (norm_by_src_rate)
    factor = _openmc_problem.tallyMultiplier(_tally_score[local_score],
                                             _local_mean_tally[local_score]) *
             _openmc_problem.scaling() * _openmc_problem.scaling();

  _modal_values.resize(_tally_score.size());
  _modal_values[local_score].assign(_num_ext_filter_bins, std::vector<Real>(n_bins, 0.0));

  for (unsigned int ext_bin = 0; ext_bin < _num_ext_filter_bins; ++ext_bin)
  {
    for (int e = 0; e < n_elem; ++e)
    {
      const int first = _surface_fet_filter->offset(e);
      const int n_modes = _surface_fet_filter->element(e).n_modes();
      if (n_modes == 0)
        continue;

      for (int i = 0; i < n_modes; ++i)
        _modal_values[local_score][ext_bin][first + i] = tally_vals[local_score](ext_bin * n_bins + first + i) * factor;

      total += _ext_bins_to_skip[ext_bin] ? 0.0 : tally_vals[local_score](ext_bin * n_bins + first);
    }
  }

  // Write the expansion into the discontinuous variable: the value at each vertex of a facet is
  // stored in that vertex's dof of the facet's element. A vertex shared by several facets of the
  // same element (an element with more than one face on the surface) gets their average.
  auto & solution = _aux.solution();
  for (unsigned int ext_bin = 0; ext_bin < _num_ext_filter_bins; ++ext_bin)
  {
    const auto var = var_numbers[_num_ext_filter_bins * local_score + ext_bin];
    std::map<dof_id_type, std::pair<Real, unsigned int>> dofs; // dof -> (sum, count)

    for (int e = 0; e < n_elem; ++e)
    {
      const auto & fs = _facet_sides[e];
      const auto * elem = _openmc_problem.getMooseMesh().queryElemPtr(fs.elem_id);
      if (!elem || !_openmc_problem.isLocalElem(elem))
        continue;

      const auto side = elem->side_ptr(fs.side);
      for (unsigned int n = 0; n < side->n_vertices(); ++n)
      {
        const unsigned int local_node = elem->get_node_index(side->node_ptr(n));
        const auto dof = elem->dof_number(_aux.number(), var, local_node);
        auto & entry = dofs[dof];
        entry.first += reconstruct(e, *side->node_ptr(n), local_score, ext_bin);
        entry.second++;
      }
    }

    for (const auto & [dof, entry] : dofs)
      solution.set(dof, entry.first / entry.second);
  }

  return total;
}

Real
FETSurfaceTally::reconstruct(int facet,
                             const libMesh::Point & p,
                             unsigned int local_score,
                             unsigned int ext_bin) const
{
  if (local_score >= _modal_values.size() || ext_bin >= _modal_values[local_score].size() ||
      _modal_values[local_score][ext_bin].empty())
    mooseError("The tally results have not been stored yet; cannot reconstruct the expansion.");
  if (facet < 0 || facet >= _surface_fet_filter->n_elements())
    mooseError("Facet index ", facet, " is out of range for the surface FET tally.");

  const auto & element = _surface_fet_filter->element(facet);
  if (element.n_modes() == 0)
    return 0.0;

  const auto r = p * _openmc_problem.scaling();

  // f(r) = sum_i N_i^2 * (tally_i) * psi_i(r) / area
  double psi[openmc::fet::MAX_MODES];
  const int n = element.eval(openmc::Position(r(0), r(1), r(2)), psi);
  const Real * values = _modal_values[local_score][ext_bin].data() + _surface_fet_filter->offset(facet);

  Real f = 0.0;
  for (int i = 0; i < n; ++i)
    f += element.mode_norm_sq(i) * values[i] * psi[i];
  return f / element.volume();
}

#endif
