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
#include "openmc/tallies/filter_sptl_legendre.h"

class LegendreTally : public FETallyBase
{
public:

    static InputParameters validParams();

    LegendreTally(const InputParameters & parameters);

    virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override;

protected:

    virtual FunctionSeries* getFunctionSeries(std::string name) override;

    virtual int getNumOrders()override {return 3;};

    virtual Real getVolume() override;

    // lower left corner of the bounding box
    std::vector<Real> _min;

    // upper right corner of the bounding box
    std::vector<Real> _max;
};
