# Aircraft Wing Spar Optimization

Mass optimization of a high performance aircraft wing spar while under loading uncertainty.

## Requirements

- MATLAB
- Optimization Toolbox (for `fmincon`)

## Quick Start

```matlab
run_opt
```

This runs the full optimization workflow:
- Visualizes the initial (nominal) spar design
- Optimizes spar geometry to minimize mass while satisfying stress constraints
- Shows convergence plots (optimality, feasibility, objective)
- Runs convergence studies to show validity of solution
- Displays final optimized geometry and stress distribution

## Key Files

- `run_opt.m` - Main script: runs optimization and generates plots
- `objective.m` - Spar mass calculation (with complex-step differentiation)
- `noncon.m` - Stress constraints (mean + 6σ < ultimate strength)
- `spar_stress_analysis.m` - Uncertainty quantification using Gauss-Hermite quadrature

## Core FEM Functions

- `CalcBeamDisplacement.m` - Solves for nodal displacements
- `CalcBeamStress.m` - Computes stresses from displacements
- `CalcElemStiff.m` - Element stiffness matrix
- `CalcElemLoad.m` - Element load vector
- `HermiteBasis.m`, `DHermiteBasis.m`, `D2HermiteBasis.m` - Cubic Hermite shape functions
- `GaussQuad.m` - Gauss quadrature integration

## Model Parameters

Default values in `run_opt.m`:

- Spar length: 7.5 m
- Material density: 1600 kg/m³
- Young's modulus: 70 GPa
- Aircraft mass: 500 kg
- Ultimate strength: 600 MPa

## Notes

- Uses Euler-Bernoulli beam theory with cubic Hermite elements
- Stress constraints use a 6σ reliability criterion
- Hollow circular cross-sections with inner/outer radii as design variables
- Complex-step differentiation for accurate gradients
