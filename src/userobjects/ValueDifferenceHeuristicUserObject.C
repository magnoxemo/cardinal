#include "ValueDifferenceHeuristicUserObject.h"
#include "MooseUtils.h"

registerMooseObject("CardinalApp", ValueDifferenceHeuristicUserObject);

InputParameters
ValueDifferenceHeuristicUserObject::validParams()
{

  InputParameters params = ClusteringUserObject::validParams();
  params.addRequiredParam<Real>("tolerance", "Maximum allowed percentage of difference"
                                             " in scores for elements to be clustered");
  params.addClassDescription("Clusters elements whose score differences"
                             " are less than the specified tolerance.");

  return params;
}

ValueDifferenceHeuristicUserObject::ValueDifferenceHeuristicUserObject(
    const InputParameters & params)
  : ClusteringUserObject(params),
    _tolerance(getParam<Real>("tolerance"))
{
}

bool
ValueDifferenceHeuristicUserObject::belongsToCluster(libMesh::Elem * base_element,
                                                     libMesh::Elem * neighbor_element)
{
  Real base_score = getMetricData(base_element);
  Real neighbor_score = getMetricData(neighbor_element);
  if (MooseUtils::isZero(base_score))
      return true;
  return std::abs((base_score-neighbor_score)/base_score) < _tolerance;
}