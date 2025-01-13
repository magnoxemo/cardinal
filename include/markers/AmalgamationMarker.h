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

#pragma once

#include "QuadraturePointMarker.h"

/**
 * A marker which takes a given variable and checks if it is within a certain range of the extrme value.
 */
class AmalgamationMarker : public QuadraturePointMarker
{
public:
  static InputParameters validParams();

  AmalgamationMarker (const InputParameters & parameters);

protected:
   virtual Marker::MarkerValue computeQpMarker() override;

  /* user defined tolerance to compare the deviation of the 
    _u[_qp] from the local average value 
     */
  const Real _tolerance;

};