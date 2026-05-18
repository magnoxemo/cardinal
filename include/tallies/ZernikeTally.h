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

#include "FETallyBase.h"
#include "openmc/tallies/filter_zernike.h"
#include <cmath>

class ZernikeTally : public FETallyBase
{
public:

    static InputParameters validParams();

    ZernikeTally(const InputParameters & parameters);

    virtual std::pair<long unsigned int, std::vector<openmc::Filter *>> 
    spatialFilters() override;

protected:

    virtual FunctionSeries* getFunctionSeries(std::string name) override;

    virtual int getNumOrders() override {return 2;};

    virtual Real getVolume() override {return 2 * _range * libMesh::pi * std::pow(_radius, 2);};

    Real _radius;

    Point _centroid;

    Real _range;
};

