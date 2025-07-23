#include "ValueFractionHeuristicUserObject.h"

registerMooseObject("CardinalApp", ValueFractionHeuristicUserObject);

InputParameters
ValueFractionHeuristicUserObject::validParams()
{

  InputParameters param = ClusteringHeuristicUserObjectBase::validParams();
  param.addRequiredParam<Real>("upper_fraction", "upper percentage of error for the heuristics");
  param.addRequiredParam<Real>("lower_fraction", "lower percentage of error for the heuristics");
  param.addClassDescription(" mimics the error fraction marker in moose");

  return param;
}
ValueFractionHeuristicUserObject::ValueFractionHeuristicUserObject(const InputParameters & params)
  : ClusteringHeuristicUserObjectBase(params),
    _upper_fraction(getParam<Real>("upper_fraction")),
    _lower_fraction(getParam<Real>("lower_fraction")),
    _max(std::numeric_limits<Real>::max()),
    _min(std::numeric_limits<Real>::min())
{
}

void
ValueFractionHeuristicUserObject::extremesFinder()
{
  double score;
  for (auto & elem : _mesh.active_element_ptr_range())
  {
    score = getMetricData(elem);
    if (_max < score)
    {
      _max = score;
    }
    if (_min > score)
    {
      _min = score;
    }
  }
  _upper_cut_off = (1 - _upper_fraction) * _max;
  _lower_cut_off = _lower_fraction * (_max - _min) + _min;
}
bool
ValueFractionHeuristicUserObject::evaluate(libMesh::Elem * base_element,
                                           libMesh::Elem * neighbor_element) const
{

  Real base_score = getMetricData(base_element);
  Real neighbor_score = getMetricData(neighbor_element);
  return (base_score > _upper_cut_off and neighbor_score > _upper_cut_off) ||
         (base_score < _lower_cut_off and neighbor_score < _lower_cut_off);
}