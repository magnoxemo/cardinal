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

#include "AmalgamationMarker.h"

registerMooseObject("MooseApp", AmalgamationMarker);

InputParameters
AmalgamationMarker ::validParams()
{
  auto params = QuadraturePointMarker::validParams();
  params.addClassDescription(
    "Mesh tally amalgamation marker based on deviation from the local mean value");
  params.addRequiredParam<Real>(
    "tolerance",
    "user defined tolerance to compare the deviation of the "
    "_u[_qp] from the local average value"
    "if (abs((_local_avg-_u[_qp])/_local_avg)<tolerance) then AMALGAMATE");
    //To do
    /*introduce _element_age concept*/
  return params;
}
AmalgamationMarker::AmalgamationMarker(const InputParameters& parameters)
  : QuadraturePointMarker(parameters), _tolerance(parameters.get<Real>("tolerance"))
{
}

Marker::MarkerValue
AmalgamationMarker::computeQpMarker()
{ 
  /* calculate local average value.
  Loop over the local elements and store the solution in a Real var*/
  Real local_avg=0;
  Real sum=0;
  for (unsigned int qp=0;qp<_qrule->n_points();qp++)
    {
        sum+=_u_neighbor[qp];
    }
  // qurle->n_points() should never be zero. But do I need to use if statement to check? 
  local_avg=sum/_qrule->n_points();
  //local_avg it could be zero. So ig it needs to be checked 
  if (std::abs(local_avg)<TOLERANCE)
  {
    return MarkerValue::DO_NOTHING;
  }
  // Avoid a divide by zero in the comparison metric and just return a do nothing state.
  if ((std::abs((local_avg-_u[_qp])/local_avg)<_tolerance))
    return MarkerValue::AMALGAMATE;
    //solution-> add the AMALGAMATE enum in moose src
  else
    return MarkerValue::DO_NOTHING;

}