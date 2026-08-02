# Master-s-code

## MATLAB Riblet Validation Tools

### `chapter4_riblet_validation_missile.m`

Processes the AGM-158 STL geometry to identify the vehicle’s principal axis and extract a representative mid body radius profile. The script creates Mach indexed wall data templates, imports CFD wall quantities, calculates parameters such as wall shear stress, skin friction coefficient, friction velocity, kinematic viscosity, and riblet spacing in wall units, and generates individual and combined plots. Its primary role is to connect the AGM-158 CAD geometry with the aerodynamic data required for initial riblet sizing and validation.

### `chapter4_riblet_validation_bomb_v8.m`

Performs the corresponding geometry and riblet sizing analysis for the GBU-57 over a Mach range of 0.80 to 1.40 in increments of 0.05. It generates STL-aligned wall data templates, estimates missing pressure and wall shear quantities when necessary, and evaluates the dimensionless riblet spacing and height parameters, `s+` and `h+`. The script also exports geometry data, run condition reports, per Mach figures, and Mach sweep summaries used to evaluate riblet sizing trends across the selected operating range. Together, these MATLAB files provide the numerical link between the STL-based CAD models, boundary layer scaling relations, and the later patterned-geometry and CFD portions of the project.
