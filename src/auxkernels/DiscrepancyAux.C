#include "DiscrepancyAux.h"
#include "SolutionUserObjectBase.h"
#include "AuxiliarySystem.h"

registerMooseObject("CardinalApp", DiscrepancyAux);

InputParameters
DiscrepancyAux::validParams(){
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("calculates the relative difference between two solution field.");
  params.addRequiredParam<UserObjectName>("solution", "The name of the SolutionUserObject");
  params.addRequiredParam<std::string>("out_mesh_variable","The name of the variable to extract from the file");
  params.addParam<AuxVariableName>("in_mesh_variable","comes from the mesh");
  return params;
}

DiscrepancyAux::DiscrepancyAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _solution_object(getUserObject<SolutionUserObjectBase>("solution")),
    _base_variable_name(getParam<AuxVariableName>("in_mesh_variable")),
    _base_variable(_subproblem.getVariable(_tid, _base_variable_name)),
    _dof_map(_aux_sys.dofMap()),
    _base_variable_index(_aux_sys.getVariable(_tid, _base_variable_name).number()),
    _from_variable_name(getParam<std::string>("out_mesh_variable"))

{
}



Real
DiscrepancyAux::computeValue()
{
  //get the score from the amr mesh
  //amr_score = get the solution from the mesh
  // The value to output
  Real base_variable_value = getBaseVariableValueByElement(_current_elem);
  Real mesh_amalgamation_output = _solution_object.pointValue(_t, _current_elem->vertex_average(),_from_variable_name );
  /*our main investigation will be how that point value get projected*/

  return (base_variable_value - mesh_amalgamation_output)/base_variable_value ;
}

Real
DiscrepancyAux::getBaseVariableValueByElement(const libMesh::Elem * elem) const
{
  std::vector<libMesh::dof_id_type> dof_indices;
  std::vector<double> solution_value(1);
  _dof_map.dof_indices(elem, dof_indices, _base_variable_index);
  _aux_sys.solution().get(dof_indices, solution_value);

  return solution_value[0];
}