#pragma once

#include "TallyBase.h"
#include "OpenMCCellAverageProblem.h"

#include "openmc/tallies/filter_surface_fet.h"

/**
 * A surface tally which expands the tally on each facet of one or more XDG surface-mesh surfaces
 * in an orthogonal functional expansion basis (OpenMC's SurfaceFETFilter). The filter has one bin
 * per (facet, mode) pair; the modal values are cached by storeResultsInner() and the expansion
 * can be evaluated at arbitrary points on a facet with reconstruct().
 *
 * The expansion is evaluated at the vertices of each facet and written into a discontinuous
 * (L2_LAGRANGE) auxvariable on the MOOSE element owning the facet, using a facet-to-(element,
 * side) map built by matching facet and side centroids.
 */
class FETSurfaceTally : public TallyBase
{
public:
  static InputParameters validParams();

  FETSurfaceTally(const InputParameters & parameters);

  virtual std::pair<unsigned int, openmc::Filter *> spatialFilter() override;

  /**
   * Initialize the tally. OpenMC decides whether a tally is a surface tally in set_scores(), which
   * needs the filters to already be set, while TallyBase sets the scores first; so the scores are
   * set again here, once the surface filter is attached.
   */
  virtual void initializeTally() override;

  /**
   * Evaluate the expansion of the most recently stored tally on a facet.
   * @param[in] facet facet index in the filter (surface after surface, XDG mesh manager order)
   * @param[in] p the point on the facet, in the units of the [Mesh]
   * @param[in] local_score index into the tally's local array of scores
   * @param[in] ext_bin the external filter bin
   * @return the tally per unit area at the point
   */
  Real reconstruct(int facet,
                   const libMesh::Point & p,
                   unsigned int local_score = 0,
                   unsigned int ext_bin = 0) const;

  /// Expansion order used on every facet.
  int order() const { return _order; }

  /// The OpenMC surface FET filter
  const openmc::SurfaceFETFilter * filter() const { return _surface_fet_filter; }

  /// Modal values, indexed by [local score][ext bin][filter bin]
  const std::vector<std::vector<std::vector<Real>>> & modalValues() const { return _modal_values; }

  /// The expansion is written to discontinuous Lagrange variables (one value per element vertex)
  virtual std::string auxVariableFamily() const override { return "L2_LAGRANGE"; }
  virtual std::string auxVariableOrder() const override { return "FIRST"; }

  /// A MOOSE element side which an OpenMC facet corresponds to
  struct FacetSide
  {
    dof_id_type elem_id;
    unsigned short side;
  };

  /// MOOSE element and side of each facet in the filter
  const std::vector<FacetSide> & facetSides() const { return _facet_sides; }

protected:
  /**
   * Build the map from the facets of the filter to MOOSE (element, side) pairs. XDG reads the mesh
   * through its own libMesh mesh, which may renumber elements and nodes, so each facet is matched
   * by geometry: its centroid (in cm) must be the centroid of a side of the [Mesh] times the
   * scaling. This requires the [Mesh] to have the same geometry as the mesh OpenMC loaded.
   */
  void buildFacetSideMap(const std::vector<int32_t> & surface_indices);

  virtual Real storeResultsInner(const std::vector<unsigned int> & var_numbers,
                                 unsigned int local_score,
                                 const std::vector<OMCTensor> & tally_vals,
                                 bool norm_by_src_rate = true) override;

  /// OpenMC IDs of the surfaces to tally on
  const std::vector<int32_t> _surface_ids;

  /// Expansion order on every facet
  const int _order;

  /// The OpenMC surface FET filter
  openmc::SurfaceFETFilter * _surface_fet_filter = nullptr;

  /// Modal values, indexed by [local score][ext bin][filter bin]
  std::vector<std::vector<std::vector<Real>>> _modal_values;

  /// MOOSE element and side of each facet, in the facet order of the filter
  std::vector<FacetSide> _facet_sides;
};
