#pragma once

#include "ClusteringHeuristicUserObjectBase.h"
class ValueRangeHeuristicUserObject : public ClusteringHeuristicUserObjectBase
{

public:
  static InputParameters validParams();
  ValueRangeHeuristicUserObject(const InputParameters & params);

protected:
  ///method for evaluting if element should be clusterd or not
  virtual bool evaluate(libMesh::Elem * base_element,
                                libMesh::Elem * neighbor_elem) const override;

  const Real _tolerance_percentage;
  Real _value;
  const Real _upper_limit;
  const Real _lower_limit;

private:
  ///determines if the score is with the range (_lower_limit, _upper_limit)
  bool isInsideTheRange(libMesh::Elem * element) const ;
};
