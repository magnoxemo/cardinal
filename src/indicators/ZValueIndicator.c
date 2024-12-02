#include "ZValueIndicator.h"

registerMooseObject("CardinalApp",ZValueIndicator);

InputParameters
ZValueIndicator::validParams()
{
    InputParameters params = ElementIntegralIndicator::validParams();
    params.addClassDescription(
        "An indicator which uses the z value of solution distribution as an error estimate"
    )
    /*params.addClassDescription(
    "An indicator which uses the integral of the norm of a vector value over an element as an error estimate.");*/
    params.addRequiredParam<VariableName>("variable",
                                        "The name of the vector variable that this indicator operates on.");
    std::vector<SubdomainName> everywhere(1, "ANY_BLOCK_ID");
    params.addParam<std::vector<SubdomainName>>(
        "block", everywhere, "block ID or name where the object works");
    return params;
}

ZValueIndicator::ZValueIndicator(const InputParameters & parameters)
: Indicator(parameters),
MooseVariableInterface<RealVectorValue>(this,
                                        false,
                                        "variable",
                                        Moose::VarKindType::VAR_ANY,
                                        Moose::VarFieldType::VAR_FIELD_VECTOR),
    _field_var(_subproblem.getStandardVariable(_tid, name())),
    _current_elem(_field_var.currentElem()),
    _current_elem_volume(_assembly.elemVolume()),
    _qrule(_assembly.qRule()),
    _JxW(_assembly.JxW()),
    _coord(_assembly.coordTransformation()),
    _vec_var(mooseVariableField()),
    _zcore_vec(_var.sln())
{
addMooseVariableDependency(&mooseVariableField());
}

void
ZValueIndicator::computeIndicator()
{
    Real heat_source_sum = 0;
    for (unsigned int qp = 0; qp < _u.size(); ++qp){
        heat_source_sum  += _u[_qp];
    }

    Real mean_heat_source = heat_source_sum/_u.size(); 
    Real std_heat_source=0;
    Real sum =0;

    for (unsigned int qp = 0; qp < _u.size(); ++qp){
        sum  += std::pow((_u[qp]-mean_heat_source),2)/__u.size();
    }

    std_heat_source= std::sqrt(sum);

    for (unsigned int qp = 0; qp < _u.size(); ++qp){
        _zscore_vec[qp]=std::abs(_u[qp]-mean_heat_source)/std_heat_source;
        //_field_var.setNodalValue(zscore);
    }
    _field_var[_current_elem]= _zscore_vec[_current_elem]

}



/* my questions:

All the examples are given for the quadrature points.
even if I follow that I need a create a same variable like _u which will hold the z_value accros the 
mesh. How do I do that? 

one way: 

1.   /// Holds the solution gradient at the current quadrature points 
const VariableValue & _u;
ref: ElementIndicator.h

or 

https://github.com/idaholab/moose/blob/next/framework/include/indicators/ElementIndicator.h




or 

_u(_var.sln()),
_grad_u(_var.gradSln())    _u(_var.sln()),
_grad_u(_var.gradSln())

https://github.com/idaholab/moose/blob/next/framework/src/indicators/ElementIndicator.C


also I am bit confused about that if the computeindicator method calculates indicator for just one _qp? 


void
VectorMagnitudeIndicator::computeIndicator()
{
Real sum = 0;
for (unsigned int qp = 0; qp < _qrule->n_points(); ++qp)
sum += _JxW[qp] * _coord[qp] * _u_vec[qp].norm();

}

https://github.com/nuclearkevin/cardinal/blob/amr_fix/src/indicators/VectorMagnitudeIndicator.C


*/