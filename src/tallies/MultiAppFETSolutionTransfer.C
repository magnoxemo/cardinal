#include "MultiAppFETSolutionTransfer.h"

#include "FETMeshTally.h"
#include "FETSurfaceTally.h"
#include "OpenMCCellAverageProblem.h"
#include "MooseMesh.h"

#include <algorithm>

registerMooseObject("CardinalApp", MultiAppFETSolutionTransfer);

InputParameters
MultiAppFETSolutionTransfer::validParams()
{
  InputParameters params = MultiAppGeneralFieldTransfer::validParams();
  params.addClassDescription(
      "Transfers a functional expansion tally by evaluating the expansion at the target points");
  params.addRequiredParam<std::string>(
      "fet_tally_object", "FET tally (FETMeshTally or FETSurfaceTally) to evaluate");
  params.addRequiredParam<std::string>("score", "Score of the tally to transfer, e.g. kappa_fission");
  params.addParam<unsigned int>(
      "ext_bin", 0, "External filter bin; only needed if the tally has external filters");

  params.addRangeCheckedParam<Real>(
      "surface_distance_tolerance",
      libMesh::TOLERANCE,
      "surface_distance_tolerance>=0",
      "Surface tallies only: target points within this distance of a tallied facet are projected "
      "onto the closest facet. Increase it for target meshes which only approximate the surface, "
      "e.g. a smooth sphere around a faceted one.");

  // The source is a tally, not a variable
  params.suppressParameter<std::vector<VariableName>>("source_variable");
  params.suppressParameter<std::vector<unsigned int>>("source_variable_components");
  params.suppressParameter<std::vector<BoundaryName>>("from_boundaries");
  params.suppressParameter<bool>("from_app_must_contain_point");
  return params;
}

MultiAppFETSolutionTransfer::MultiAppFETSolutionTransfer(const InputParameters & parameters) noexcept(false)
  : MultiAppGeneralFieldTransfer(parameters),
    _tally_name(getParam<std::string>("fet_tally_object")),
    _score(getParam<std::string>("score")),
    _ext_bin(getParam<unsigned int>("ext_bin")),
    _surface_tol(getParam<Real>("surface_distance_tolerance"))
{
  if (_to_var_names.size() != 1)
    paramError("variable", "Exactly one target variable is supported");
}

std::string
MultiAppFETSolutionTransfer::getDataSourceName(unsigned int var_index) const
{
  return "FET tally '" + _tally_name + "'";
}

void
MultiAppFETSolutionTransfer::prepareEvaluationOfInterpValues(const unsigned int var_index)
{
  _local_bboxes.clear();
  if (_use_bounding_boxes)
    extractLocalFromBoundingBoxes(_local_bboxes);

  _point_locators.clear();
  _mesh_tallies.assign(_from_problems.size(), nullptr);
  _surface_tallies.assign(_from_problems.size(), nullptr);

  for (const auto i : index_range(_from_problems))
  {
    auto * problem = dynamic_cast<OpenMCCellAverageProblem *>(_from_problems[i]);
    if (!problem)
      mooseError("The source problem must be an OpenMCCellAverageProblem");

    const auto * tally = problem->getTally(_tally_name);
    if (!tally)
      paramError("fet_tally_object", "Tally '", _tally_name, "' does not exist");

    // OpenMC spells scores with dashes; accept the underscore spelling used in the input files
    std::string score = _score;
    std::replace(score.begin(), score.end(), '_', '-');
    const auto & scores = tally->getScores();
    const auto it = std::find(scores.begin(), scores.end(), score);
    if (it == scores.end())
      paramError("score", "Tally '", _tally_name, "' does not have score '", _score, "'");
    _local_score = std::distance(scores.begin(), it);

    _mesh_tallies[i] = dynamic_cast<const FETMeshTally *>(tally);
    _surface_tallies[i] = dynamic_cast<const FETSurfaceTally *>(tally);
    if (!_mesh_tallies[i] && !_surface_tallies[i])
      paramError("fet_tally_object", "Tally '", _tally_name, "' is not a FET tally");

    auto pl = _from_problems[i]->mesh().getMesh().sub_point_locator();
    pl->enable_out_of_mesh_mode();
    _point_locators.push_back(std::move(pl));
  }
}

namespace
{
/// Closest point to p on the segment (a, b)
Point
closestOnSegment(const Point & p, const Point & a, const Point & b)
{
  const Point ab = b - a;
  const Real t = std::max(0.0, std::min(1.0, ((p - a) * ab) / ab.norm_sq()));
  return a + ab * t;
}

/// Closest point to p on the triangle (a, b, c)
Point
closestOnTriangle(const Point & p, const Point & a, const Point & b, const Point & c)
{
  const Point n = (b - a).cross(c - a);
  const Point q = p - n * (((p - a) * n) / n.norm_sq());

  // the projection onto the plane is the answer if it lies inside all three edges
  if ((b - a).cross(q - a) * n >= 0 && (c - b).cross(q - b) * n >= 0 &&
      (a - c).cross(q - c) * n >= 0)
    return q;

  Point best = closestOnSegment(p, a, b);
  for (const auto & candidate : {closestOnSegment(p, b, c), closestOnSegment(p, c, a)})
    if ((candidate - p).norm_sq() < (best - p).norm_sq())
      best = candidate;
  return best;
}
}

int
MultiAppFETSolutionTransfer::closestFacet(const unsigned int i_from,
                                          const Point & p,
                                          Point & p_facet) const
{
  const auto & mesh = _from_problems[i_from]->mesh();
  const auto & fs = _surface_tallies[i_from]->facetSides();

  int facet = -1;
  Real best = _surface_tol;
  for (const auto f : index_range(fs))
  {
    const auto side = mesh.getMesh().elem_ptr(fs[f].elem_id)->side_ptr(fs[f].side);
    // triangles, and (planar) quadrilaterals split into the triangles (0, 1, 2) and (0, 2, 3)
    for (unsigned int v = 2; v < side->n_vertices(); ++v)
    {
      const auto q =
          closestOnTriangle(p, side->point(0), side->point(v - 1), side->point(v));
      const Real distance = (q - p).norm();
      if (distance <= best)
      {
        best = distance;
        facet = f;
        p_facet = q;
      }
    }
  }
  return facet;
}

Real
MultiAppFETSolutionTransfer::evaluateTally(const unsigned int i_from, const Point & p) const
{
  if (_mesh_tallies[i_from])
  {
    const auto * elem = (*_point_locators[i_from])(p);
    if (!elem)
      return GeneralFieldTransfer::OutOfMeshValue;
    return _mesh_tallies[i_from]->reconstruct(elem, p, _local_score, _ext_bin);
  }

  Point p_facet;
  const int facet = closestFacet(i_from, p, p_facet);
  if (facet < 0)
    return GeneralFieldTransfer::OutOfMeshValue;

  return _surface_tallies[i_from]->reconstruct(facet, p_facet, _local_score, _ext_bin);
}

void
MultiAppFETSolutionTransfer::evaluateInterpValues(
    const unsigned int,
    const std::vector<std::pair<Point, unsigned int>> & incoming_points,
    std::vector<std::pair<Real, Real>> & outgoing_vals)
{
  dof_id_type i_pt = 0;
  for (auto & [pt, mesh_div] : incoming_points)
  {
    outgoing_vals[i_pt] = {GeneralFieldTransfer::OutOfMeshValue,
                           GeneralFieldTransfer::OutOfMeshValue};

    for (const auto i_from : index_range(_from_problems))
    {
      Real distance = 0;
      if (!acceptPointInOriginMesh(i_from, _local_bboxes, pt, mesh_div, distance))
        continue;

      const auto local_pt = getPointInSourceAppFrame(pt, i_from, "FET evaluation");
      const Real val = evaluateTally(i_from, local_pt);
      if (val == GeneralFieldTransfer::OutOfMeshValue)
        continue;

      if (distance < outgoing_vals[i_pt].second)
        outgoing_vals[i_pt] = {val, distance};
    }
    ++i_pt;
  }
}
