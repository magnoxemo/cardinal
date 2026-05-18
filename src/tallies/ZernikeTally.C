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
#include "ZernikeTally.h"
#include "openmc/tallies/filter_zernike.h"

registerMooseObject("CardinalApp", ZernikeTally);

InputParameters
ZernikeTally::validParams()
{
  auto params = FETallyBase::validParams();
  params.addClassDescription("A class which implements Zernike expansion tallies.");
  params.addRequiredParam<Point>("centroid",
                                 "The centroid of the expansion bounds (x, y, z) "
                                 "for the bounding box.");
  params.addRequiredParam<Real>("radius",
                                "The radius of the Zernike polynomials.");
  params.addRequiredParam<Real>("half_range",
                                "Distance from centroid to upper/ lower bound.");
  return params;
}

ZernikeTally::ZernikeTally(const InputParameters & parameters)
  : FETallyBase(parameters),
    _radius(getParam<Real>("radius")),
    _centroid(getParam<Point>("centroid")),
    _range(getParam<Real>("half_range"))
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
ZernikeTally::spatialFilters()
{
  std::vector<openmc::Filter*> filters;
  //Legendre
  auto * lfilter = dynamic_cast<openmc::SpatialLegendreFilter *>(openmc::Filter::create("spatiallegendre"));
  lfilter->set_minmax(_centroid(2) - _range, _centroid(2) + _range);
  lfilter->set_order(_orders.at(1));
  lfilter->set_axis(openmc::LegendreAxis::z);
  filters.push_back(static_cast<openmc::Filter *>(lfilter));
  //Zernike
  auto zfilter = dynamic_cast<openmc::ZernikeFilter *>(openmc::Filter::create("zernike"));
  zfilter->set_order(_orders.at(0));
  zfilter->set_x(_centroid(0));
  zfilter->set_y(_centroid(1));
  zfilter->set_r(_radius);
  filters.push_back(static_cast<openmc::Filter *>(zfilter));

  return std::make_pair(openmc::model::tally_filters.size() - 2, filters);
}

FunctionSeries*
ZernikeTally::getFunctionSeries(std::string name)
{
  std::vector<Real> bounds{_centroid(2)-_range, _centroid(2)+_range,
                           _centroid(0), _centroid(1), _radius};
  std::vector<unsigned> orders{_orders.at(1), _orders.at(0)};

  return _openmc_problem.makeFunctionSeries(name, "CylindricalDuo", orders, bounds);
}

#endif