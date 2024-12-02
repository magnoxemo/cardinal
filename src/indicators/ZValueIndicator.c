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
