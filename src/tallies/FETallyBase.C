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

#include "FETallyBase.h"

InputParameters
FETallyBase::validParams()
{
  auto params = TallyBase::validParams();
  params.addRequiredParam<std::vector<unsigned>>("orders",
                                                 "The orders for the expansions in each dimension. "
                                                 "These must be (x, y, z) for Legendre and "
                                                 "(rtheta, z) for Zernike.");
  params.addParam<std::string>("function_suffix","_function",
                               "The suffix to append to the score as the name of the"
                               "function object holding the Legendre expansion.");
  params.set<MultiMooseEnum>("output") = "UNRELAXED_TALLY";
  return params;
}

FETallyBase::FETallyBase(const InputParameters & parameters)
  : TallyBase(parameters),
    _orders(getParam<std::vector<unsigned>>("orders")),
    _function_suffix(getParam<std::string>("function_suffix"))
{

  //overriding auxvariable names, don't want to create any
  if (isParamValid("name"))
  {
    mooseWarning(this->_name + " does not have any ElementalAuxVariables "
                               "associated with it! "+this->_name+" creates functions, "
                               "the names of which are controllable by \"function_suffix\". "
                               "Clearing \"name\" parameter...");
  }
  _tally_name.clear();

  // initializing number of functions and coefficients
  _functions.resize(_tally_score.size());

  /**
     * OpenMC spatial FETs only support the collision estimator
   */
  if (isParamValid("estimator"))
  {
    if (_estimator != openmc::TallyEstimator::COLLISION)
      paramError("estimator",
                 "Collision estimators are currently the only compatible "
                 "estimator type for Spatial Legendre expansion tallies!");
  }
  else
    _estimator = openmc::TallyEstimator::COLLISION;

  // initializing the first moments
  _first_moments.resize(_tally_score.size());

}


void
FETallyBase::resetTally()
{
  // Erase the tally.
  openmc::model::tallies.erase(openmc::model::tallies.begin() + _local_tally_index);

  // Erase the global normalization tally.
  if (addingGlobalTally())
    openmc::model::tallies.erase(openmc::model::tallies.begin() + _global_tally_index);

  for (int i = getNumOrders() - 1; i >= 0; --i)
    openmc::model::tally_filters.erase(openmc::model::tally_filters.begin() + _filter_index + i);
}

void
FETallyBase::computeSumAndMean()
{

  for (unsigned int score = 0; score < _tally_score.size(); ++score)
  {
    auto openmc_coeffs = _openmc_problem.tallySum(_local_tally, score);

    _first_moments[score] = openmc_coeffs(0);

    _local_sum_tally[score] = _first_moments[score];
    _local_mean_tally[score] = _first_moments[score] / getVolume();

    if (addingGlobalTally())
      _global_sum_tally[score] = _openmc_problem.tallySumAcrossBins({_global_tally}, score);
  }
}

Real FETallyBase::storeResultsInner(const std::vector<unsigned int> & var_numbers,
                               unsigned int local_score,
                               const std::vector<OMCTensor> & tally_vals,
                               bool norm_by_src_rate)
{
  // tally_vals[local_score] contains the normalized FET coefficients (divided by the
  // normalization factor in relaxAndNormalizeTally). Copy them into the FunctionSeries
  // so the expansion is available as a MOOSE Function object.
  std::vector<Real> coeffs;
  coeffs.reserve(tally_vals[local_score].size());
  for (std::size_t i = 0; i < tally_vals[local_score].size(); ++i)
    coeffs.push_back(tally_vals[local_score](i));

  _functions[local_score]->setCoefficients(coeffs);

  // The zeroth-order coefficient of the normalized expansion equals the
  // normalized integral over the domain, which must be ≈ 1.0 for a properly
  // constructed FET tally. Return it so that TallyBase::checkNormalization
  // can verify the expansion sums correctly.
  return coeffs.empty() ? 0.0 : coeffs[0];
}

#endif