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

#include "TallyBase.h"
#include "OpenMCCellAverageProblem.h"
#include "FunctionSeries.h"

#include "openmc/tensor.h"

class FETallyBase : public TallyBase
{
public:

  static InputParameters validParams();

  FETallyBase(const InputParameters & parameters);


  virtual void resetTally() override;

  virtual void computeSumAndMean() override;

  virtual std::pair<unsigned int, openmc::Filter *> spatialFilter()=0;

protected:
  virtual Real getNormalizedCoefficients(std::size_t index) const = 0;

  virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 const std::vector<OMCTensor> & tally_vals,
                                 bool norm_by_src_rate = true) override;

  virtual FunctionSeries* getFunctionSeries(std::string name) = 0;

  virtual int getNumOrders(){return 0;};
  virtual Real getVolume() = 0;

  std::vector<unsigned int> _orders;

  std::vector<FunctionSeries*> _functions;
  std::vector<Real> _first_moments;

  std::string _function_suffix;

};