#include "DiscrepenyAux.h"
#include "MooseError.h"
#include "SolutionUserObjectBase.h"

registerMooseObject("CardinalApp", DiscrepancyAux);

InputParameters
DiscrepancyAux::validParams(){
  InputParameters params = AuxKernel::validParams();
  params.addRequiredParam<VariableName>("base_variable_name","name of the variable which's relative"
                                                              " difference will be calculated");
  params.addRequiredParam<UserObjectName>("solution", "The name of the SolutionUserObject");
  params.addParam<std::string>("from_variable","The name of the variable to extract from the file");

  params.addClassDescription("calculates the relative difference between two solution field.");

  return params;
}

DiscrepancyAux::DiscrepancyAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    MooseVariableInterFace<Real>(this, false, "base_variable_name")
    _solution_object(getUserObject<SolutionUserObjectBase>("solution")),

{
}

void
DiscrepancyAux::initialSetup()
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
