#pragma once

#include "ClusteringHeuristicUserObjectBase.h"

class ValueDifferenceHeuristicUserObject : public ClusteringHeuristicUserObjectBase
{

public:
  static InputParameters validParams();
  ValueDifferenceHeuristicUserObject(const InputParameters & params);

protected:
  /// method for evaluting if element should be clusterd or not
  virtual bool evaluate(libMesh::Elem * base_element,
                                libMesh::Elem * neighbor_elem)const override;

  /// relative tolerance
  const Real _tolerance;
};