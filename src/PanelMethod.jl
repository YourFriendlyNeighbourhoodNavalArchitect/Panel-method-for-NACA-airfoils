# Constant-strength source/vortex panel method for 2D airfoils in incompressible potential flow.
# Panels run counter-clockwise. Therefore, γ is positive clockwise.

module PanelMethod

using LinearAlgebra: LinearAlgebra

export readAirfoil, naca4,
       Panels, discretizeAirfoil, nPanels, chordLength, upperSurface, lowerSurface,
       sourceVelocity, influenceMatrices,
       PanelSolution, solvePanels,
       liftCoefficient, pressureLiftCoefficient, momentCoefficient,
       flowVelocity, velocityField, isInside

include("geometry.jl")
include("influence.jl")
include("solver.jl")
include("postprocess.jl")

end
