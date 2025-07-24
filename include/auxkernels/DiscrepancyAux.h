#pragma once
#include "AuxKernel.h"

//forward declaration
class SolutionUserObjectBase;

class DiscrepancyAux : public AuxKernel
{
public:
  static InputParameters validParams();

  DiscrepancyAux(const InputParameters & parameters);
  virtual void initialSetup() override{};

protected:

  virtual Real computeValue() override;

  /// Reference to the SolutionUserObject storing the solution
  const SolutionUserObjectBase & _solution_object;

  Real getBaseVariableValueByElement(const libMesh::Elem* element) const;

  /// The variable name of interest
  std::string _mesh_amalgamation_var_name;

  const AuxVariableName & _base_variable_name;
  MooseVariableBase &  _base_variable;

  ///DOF map
  libMesh::DofMap & _dof_map;

  ///Metric variable index
  const unsigned int _base_variable_index;
  std::string _from_variable_name;

};