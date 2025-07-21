#pragma once
#include "AuxKernel.h"
//forward declaration
class SolutionUserObjectBase;

class DiscrepencyAux : public AuxKernel
{
public:
  static InputParameters validParams();

  DiscrepencyAux(const InputParameters & parameters);
  virtual void initialSetup() override;

protected:

  virtual Real computeValue() override;

  /// Reference to the SolutionUserObject storing the solution
  const SolutionUserObjectBase & _solution_object;

  Real getBaseVariableValueByElement(libMesh::Elem* element) const;

  /// The variable name of interest
  std::string _mesh_amalgamation_var_name;

  AuxVariableName & _base_variable_name;
  MooseVariableBase &  _base_variable;

  ///AuxiliarySystem reference
  AuxiliarySystem & _auxiliary_system;

  ///DOF map
  libMesh::DofMap & _dof_map;

  ///Metric variable index
  const unsigned int _base_variable_index;

};