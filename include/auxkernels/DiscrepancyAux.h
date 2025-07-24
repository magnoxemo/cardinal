#pragma once
#include "AuxKernel.h"
#include "MooseVariableInterFace.h"
//forward declaration
class SolutionUserObjectBase;

class DiscrepencyAux : public AuxKernel, public MooseVariableInterface<Real>
{
public:
  static InputParameters validParams();

  DiscrepencyAux(const InputParameters & parameters);
  virtual void initialSetup() override;

protected:

  virtual Real computeValue() override;

  /// Reference to the SolutionUserObject storing the solution
  const SolutionUserObjectBase & _solution_object;

  //base variable value
  const VariableValue & _base_variable_value;


  /// The variable name of interest
  std::string _mesh_amalgamation_var_name;


};