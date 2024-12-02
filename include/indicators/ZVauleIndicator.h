#pragma once

#include "Indicator.h"
#include "MooseVariableInterface.h"

class ZValueIndicator : public ElementIntegralIndicator, 
                        public MooseVariableInterface<RealVectorValue>,
                        public MooseVariableInterface<Real>
{
public:
    static InputParameters validParams();
    ZValueIndicator(const InputParameters & parameters);
    virtual void computeIndicator() override;

protected:
    MooseVariable & _field_var;

    const Elem * const & _current_elem;
    /// Volume of the current element
    const Real & _current_elem_volume;

    const QBase * const & _qrule;
    const MooseArray<Real> & _coord;
    
    /// Holds the solution at current quadrature points
    MooseVariableField<RealVectorValue> & _var;

    const VectorVariableValue & _u;
    /// Holds the zvalue at current quadrature points
    const VectorVariableValue & _zscore_vec;
}