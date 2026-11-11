#pragma once

#include "MeshTally.h"

#include "openmc/tallies/filter_mesh_fet.h"

/**
 * An unstructured mesh tally which expands the tally within each (linear tetrahedral) element
 * in an orthogonal functional expansion basis. The expansion is evaluated at the
 * vertices of each element and written to discontinuous linear (L2_LAGRANGE) auxvariables, and the
 * full expansion can be evaluated at arbitrary points with reconstruct().
 */
class FETMeshTally : public MeshTally
{
public:
  static InputParameters validParams();

  FETMeshTally(const InputParameters & parameters);

  /**
   * Evaluate the expansion of the most recently stored tally at a point.
   * @param[in] elem the element containing the point
   * @param[in] p the point, in the units of the [Mesh]
   * @param[in] local_score index into the tally's local array of scores
   * @param[in] ext_bin the external filter bin
   * @return the tally per unit volume at the point
   */
  Real reconstruct(const libMesh::Elem * elem,
                   const libMesh::Point & p,
                   unsigned int local_score = 0,
                   unsigned int ext_bin = 0) const;

  /// The expansion is written to discontinuous linear Lagrange variables (one value per vertex)
  virtual std::string auxVariableFamily() const override { return _family; }
  virtual std::string auxVariableOrder() const override { return _nodal_variable_order; }

  /// Expansion order used on every element.
  int order() const { return _order; }

protected:
  /// reconstruct() needs the element-to-bin map whenever block restriction is used
  virtual bool needsElementToBinMap() const override { return _use_dof_map; }

  virtual openmc::MeshFilter * createMeshFilter() override;
  virtual void configureMeshFilter() override;

  /**
   * Stores the expansion evaluated at each element vertex in the aux variables, and caches the
   * normalized modal values for reconstruct().
   */
  virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 const std::vector<OMCTensor> & tally_vals,
                                 bool norm_by_src_rate = true) override;

  /// Expansion order on every element
  /// This doesn't need to be uniform. Ig my future work is gonna be figuring out how to
  /// not make this uniform
  const int _order;

  /// The OpenMC FET mesh filter
  /// We need that to control local order of the basis function
  openmc::MeshFETFilter * _mesh_fet_filter = nullptr;

  /// Normalized modal values, indexed by [local score][ext bin][filter bin]
  std::vector<std::vector<std::vector<Real>>> _modal_values;

  const std::string _family;
  const std::string _nodal_variable_order;
};
