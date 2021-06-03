[Problem]
  type = NekRSProblem
[]

[Mesh]
  type = NekRSMesh
  boundary = '1 2 3 4 5 6 7 8'
[]

[Executioner]
  type = Transient

  [TimeStepper]
    type = NekTimeStepper
  []
[]

[Outputs]
  [out]
    type = CSV
    hide = 'flux_integral'
    execute_on = 'final'
  []
[]

[Postprocessors]
  [max_temp_side1]
    type = NekSideExtremeValue
    field = temperature
    boundary = '1'
    value_type = max
  []
  [max_temp_side2]
    type = NekSideExtremeValue
    field = temperature
    boundary = '2'
    value_type = max
  []
  [max_temp_side3]
    type = NekSideExtremeValue
    field = temperature
    boundary = '3'
    value_type = max
  []
  [max_temp_side4]
    type = NekSideExtremeValue
    field = temperature
    boundary = '4'
    value_type = max
  []
  [max_temp_side5]
    type = NekSideExtremeValue
    field = temperature
    boundary = '5'
    value_type = max
  []
  [max_temp_side6]
    type = NekSideExtremeValue
    field = temperature
    boundary = '6'
    value_type = max
  []
  [max_temp_side7]
    type = NekSideExtremeValue
    field = temperature
    boundary = '7'
    value_type = max
  []
  [max_temp_side8]
    type = NekSideExtremeValue
    field = temperature
    boundary = '8'
    value_type = max
  []
  [min_temp_side1]
    type = NekSideExtremeValue
    field = temperature
    boundary = '1'
    value_type = min
  []
  [min_temp_side2]
    type = NekSideExtremeValue
    field = temperature
    boundary = '2'
    value_type = min
  []
  [min_temp_side3]
    type = NekSideExtremeValue
    field = temperature
    boundary = '3'
    value_type = min
  []
  [min_temp_side4]
    type = NekSideExtremeValue
    field = temperature
    boundary = '4'
    value_type = min
  []
  [min_temp_side5]
    type = NekSideExtremeValue
    field = temperature
    boundary = '5'
    value_type = min
  []
  [min_temp_side6]
    type = NekSideExtremeValue
    field = temperature
    boundary = '6'
    value_type = min
  []
  [min_temp_side7]
    type = NekSideExtremeValue
    field = temperature
    boundary = '7'
    value_type = min
  []
  [min_temp_side8]
    type = NekSideExtremeValue
    field = temperature
    boundary = '8'
    value_type = min
  []

  [max_pressure_side1]
    type = NekSideExtremeValue
    field = pressure
    boundary = '1'
  []
  [max_pressure_side2]
    type = NekSideExtremeValue
    field = pressure
    boundary = '2'
  []
  [max_pressure_side3]
    type = NekSideExtremeValue
    field = pressure
    boundary = '3'
  []
  [max_pressure_side4]
    type = NekSideExtremeValue
    field = pressure
    boundary = '4'
  []
  [max_pressure_side5]
    type = NekSideExtremeValue
    field = pressure
    boundary = '5'
  []
  [max_pressure_side6]
    type = NekSideExtremeValue
    field = pressure
    boundary = '6'
  []
  [max_pressure_side7]
    type = NekSideExtremeValue
    field = pressure
    boundary = '7'
  []
  [max_pressure_side8]
    type = NekSideExtremeValue
    field = pressure
    boundary = '8'
  []
  [min_pressure_side1]
    type = NekSideExtremeValue
    field = pressure
    boundary = '1'
    value_type = min
  []
  [min_pressure_side2]
    type = NekSideExtremeValue
    field = pressure
    boundary = '2'
    value_type = min
  []
  [min_pressure_side3]
    type = NekSideExtremeValue
    field = pressure
    boundary = '3'
    value_type = min
  []
  [min_pressure_side4]
    type = NekSideExtremeValue
    field = pressure
    boundary = '4'
    value_type = min
  []
  [min_pressure_side5]
    type = NekSideExtremeValue
    field = pressure
    boundary = '5'
    value_type = min
  []
  [min_pressure_side6]
    type = NekSideExtremeValue
    field = pressure
    boundary = '6'
    value_type = min
  []
  [min_pressure_side7]
    type = NekSideExtremeValue
    field = pressure
    boundary = '7'
    value_type = min
  []
  [min_pressure_side8]
    type = NekSideExtremeValue
    field = pressure
    boundary = '8'
    value_type = min
  []

  [flux_integral]
    type = Receiver
  []
[]