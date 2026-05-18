/********************************************************************/
/*                  SOFTWARE COPYRIGHT NOTIFICATION                 */
/*                             Cardinal                             */
/*                                                                  */
/*                  (c) 2021 UChicago Argonne, LLC                  */
/*                        ALL RIGHTS RESERVED                       */
/*                                                                  */
/*                 Prepared by UChicago Argonne, LLC                */
/*               Under Contract No. DE-AC02-06CH11357               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*             Prepared by Battelle Energy Alliance, LLC            */
/*               Under Contract No. DE-AC07-05ID14517               */
/*                With the U. S. Department of Energy               */
/*                                                                  */
/*                 See LICENSE for full restrictions                */
/********************************************************************/

#ifdef ENABLE_OPENMC_COUPLING
#include "LegendreTally.h"

#include "openmc/tallies/filter_sptl_legendre.h"

registerMooseObject("CardinalApp", LegendreTally);

InputParameters
LegendreTally::validParams()
{
  auto params = FETallyBase::validParams();
  params.addClassDescription("A class which implements Legendre expansion tallies.");
  params.addRequiredParam<Point>("minimum",
                                 "The minimum bounds (x, y, z) for the bounding box.");
  params.addRequiredParam<Point>("maximum",
                                 "The maximum bounds (x, y, z) for the bounding box.");
  return params;
}

LegendreTally::LegendreTally(const InputParameters & parameters)
  : FETallyBase(parameters),
    _min(getParam<Point>("minimum")),
    _max(getParam<Point>("maximum"))
{
  //Verifying number of orders passed is correct
  if (_orders.size() != getNumOrders())
    mooseError("Cardinal only supports 3-D, the length of \"orders\" "
               "for "+this->_name+" must be equal to "
               +std::to_string(getNumOrders())+".");

  // initializing functions
  for (std::size_t index = 0; index < _tally_score.size(); ++index)
  {
    _functions.at(index) = this->getFunctionSeries(_tally_score.at(index) + _function_suffix);
  }

}

std::pair<long unsigned int, std::vector<openmc::Filter *>>
LegendreTally::spatialFilters()
{

  std::vector<openmc::Filter*> filters;

  for (int i = 0; i < 3; ++i)
  {
    auto filter = dynamic_cast<openmc::SpatialLegendreFilter *>(openmc::Filter::create("spatiallegendre"));
    filter->set_minmax(_min(i), _max(i));
    filter->set_order(_orders.at(i));
    filter->set_axis(static_cast<openmc::LegendreAxis>(i));
    filters.push_back(static_cast<openmc::Filter *>(filter));
  }

  return std::make_pair(openmc::model::tally_filters.size() - 3, filters);
}

FunctionSeries*
LegendreTally::getFunctionSeries(std::string name)
{
  std::vector<Real> _bounds{_min(0), _max(0), _min(1), _max(1), _min(2), _max(2)};

  return _openmc_problem.makeFunctionSeries(name, "Cartesian", _orders, _bounds);
}

#endif