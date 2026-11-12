#pragma once

#include "MultiAppGeneralFieldTransfer.h"
#include "libmesh/point_locator_base.h"

class FETMeshTally;
class FETSurfaceTally;

/**
 * Transfers a functional expansion tally to the nodes of a target (sub)app mesh by evaluating the
 * expansion at the target points, instead of transferring an auxvariable. The point collection,
 * bounding box checks and parallel communication come from MultiAppGeneralFieldTransfer; only the
 * evaluation of the source is specialised.
 */
class MultiAppFETSolutionTransfer : public MultiAppGeneralFieldTransfer
{
public:
  static InputParameters validParams();

  MultiAppFETSolutionTransfer(const InputParameters & parameters);

protected:
  virtual void prepareEvaluationOfInterpValues(const unsigned int var_index) override;

  virtual void
  evaluateInterpValues(const unsigned int var_index,
                       const std::vector<std::pair<Point, unsigned int>> & incoming_points,
                       std::vector<std::pair<Real, Real>> & outgoing_vals) override;

  virtual std::string getDataSourceName(unsigned int var_index) const override;

private:
  /// Evaluate the tally of source problem i_from at a point in that app's frame
  Real evaluateTally(const unsigned int i_from, const Point & p) const;

  /// Name of the FET tally object in the OpenMCCellAverageProblem
  const std::string _tally_name;
  /// Name of the score to transfer
  const std::string _score;
  /// External filter bin (0 if the tally has no external filters)
  const unsigned int _ext_bin;
  /// Index of the score within the tally; resolved when preparing the evaluation
  unsigned int _local_score = 0;

  std::vector<libMesh::BoundingBox> _local_bboxes;
  std::vector<std::unique_ptr<libMesh::PointLocatorBase>> _point_locators;
  /// Tally per source problem (only one of each pair is set)
  std::vector<const FETMeshTally *> _mesh_tallies;
  std::vector<const FETSurfaceTally *> _surface_tallies;
};
