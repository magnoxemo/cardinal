[Mesh]
  [sphere]
    type = FileMeshGenerator
    file = ../meshes/pincell.e
  []
[]

[Problem]
    type = OpenMCCellAverageProblem
    verbose = true
    power = 500.0
    check_tally_sum = false
    normalize_by_global_tally = true

    [Tallies]
        [Zernike]
            type = ZernikeTally
            orders = '10 10'
            centroid = '0 0 5'
            radius = .4125
            half_range = 5.
        []
    []
[]


[AuxVariables]
  [kappa_fission]
  []
[]

[AuxKernels]
  [FunctionAux]
    type = FunctionAux
    variable = kappa_fission
    function = kappa-fission_function
  [../]
[]

[Executioner]
    type = Transient
    num_steps = 5
[]

[Outputs]
  exodus = true
[]