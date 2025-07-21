#include "DiscrepenyAux.h"
#include "MooseError.h"
#include "SolutionUserObjectBase.h"

registerMooseObject("DiscrepencyAux", CardinalApp);

InputParameters
DiscrepencyAux::validParams(){
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<VariableName>("base_variable_name","name of the variable which's relative"
                                                              " difference will be calculated");
  params.addClassDescription("calculates the relative difference between two solution field.");
  params.addRequiredParam<UserObjectName>("solution", "The name of the SolutionUserObject");
  params.addParam<std::string>("from_variable","The name of the variable to extract from the file");
  return params;
}

DiscrepencyAux::DiscrepencyAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _solution_object(getUserObject<SolutionUserObjectBase>("solution")),
    _base_variable_name(getParam<AuxVariableName>("metric_variable_name")),
    _base_variable(_fe_problem.getVariable(_tid, _base_variable_name)),
    _auxiliary_system(_fe_problem.getAuxiliarySystem()),
    _dof_map(_auxiliary_system.dofMap()),
    _base_variable_index(_auxiliary_system.getVariable(_tid, _base_variable_name).number())

{
}

void
DiscrepencyAux::initialSetup()
{
  // If 'from_variable' is supplied, use the value
  if (isParamValid("from_variable"))
    _var_name = getParam<std::string>("from_variable");
  else
    //throw a moose error
    mooseError("We need the mesh_amalgamation score");
}

Real
DiscrepencyAux::computeValue()
{
  //get the score from the amr mesh
  //amr_score = get the solution from the mesh
  // The value to output
  Real base_variable_value = getBaseVariableValueByElement(_current_elem);
  Real mesh_amalgamtion_output = _solution_object.pointValue(_t, _current_elem->vertex_average(), _var_name);
  /*our main investigation will be how that point value get projected*/

  return (base_variable_value - mesh_amalgamation_output)/base_variable_value ;
}

Real
DiscrepencyAux::getBaseVariableValueByElement(const libMesh::Elem * elem) const
{
  std::vector<libMesh::dof_id_type> dof_indices;
  std::vector<double> solution_value(1);
  _dof_map.dof_indices(elem, dof_indices, _metric_variable_index);
  _auxiliary_system.solution().get(dof_indices, solution_value);

  return solution_value[0];
}