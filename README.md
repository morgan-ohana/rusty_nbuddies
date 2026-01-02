# Rusty Nbuddies

The program is used by running the binary with the path to the parameters file as it's only command line argument. The parameters file should be in TOML and consist of 3 tables [Initial-Conditions] [Simulation] and [Diagnostics]. These correspond to the three largely independent functions of the codebase. Any of these tables may be ommited to signal you do not want the code to do those tasks (obviously the simulation cannot run without an inital condition though). For everything in this file mass unit is solar masses, distance is kiloparsec, and time is in gigayear.

## [Initial-Conditions]

The initial conditions table must contain a **type** key which matches one of the following: ["File", "Binary", "Plummer", "NFW", "Alpha-Beta-Gamma"]. For all except "File" and "Binary" this will call the built in eddington inverter to generate your initial conditions.

### type = "File"

This is the simplest type and used for simply directing the code to load an existing initial conditions file. In this case the only other entry you'll need in your initial conditions table is `location = path/to/file`

### type = "Binary"

This creates a binary initial condition, the required additional keys in your table are then as follows. If the default value is listed as N/A the value must be specified by the user

| Key | Default | Description |
|-----|---------|-------------|
| output-directory | N/A | where to write initial conditions file to |
| m1 | N/A | mass of first particle |
| m2 | N/A | mass of second particle |
| r_min | N/A | closest approach |
| r_init | r_min | initial seperation |
| eccentricity | 0 | eccentricity of the orbit |

### type = "Plummer"

Uses the eddington inverter to populate a plummer sphere. The additional keys are :

| Key | Default | Description |
|-----|---------|-------------|
| output-directory | N/A | where to write initial conditions file to |
| N | N/A | the number of particles |
| r_s | N/A | the scale radius of the plummer sphere |
| M | N/A | total mass of plummer sphere |

### type = "NFW"

Uses the eddington inverter to populate a truncated NFW profile which exponentially decays outside some r_cutoff following the HALOGEN4MUSE prescription. The additional keys are :

| Key | Default | Description |
|-----|---------|-------------|
| output-directory | N/A | where to write initial conditions file to |
| N | N/A | the number of particles |
| r_s | N/A | the scale radius of the NFW profile |
| rho_s | N/A | the scale density of the NFW profile |
| r_cutoff | N/A | the cutoff radius of the NFW profile |

### type = "Alpha-Beta-Gamma"

Uses the eddington inverter to populate an Alpha-Beta-Gamma profile following the prescription of HALOGEN4MUSE. If beta <= 3 the total mass would be divergent so a cutoff radius outside of which it will exponentially decay is required. The additional keys are :

| Key | Default | Description |
|-----|---------|-------------|
| output-directory | N/A | where to write initial conditions file to |
| N | N/A | the number of particles |
| alpha | N/A | the alpha parameter |
| beta | N/A | the beta parameter |
| gamma | N/A | the gamma parameter |
| r_s | N/A | the scale radius of the NFW profile |
| rho_s | N/A | the scale density of the NFW profile |
| r_cutoff | N/A (Optional if beta > 3)! | the cutoff radius of the NFW profile |

## [Simulation]

The simulation table dictates how the simulation will be run. The keys for this table are listed below, for some detailed descriptions are after the table:

| Key | Default | Description |
|-----|---------|-------------|
| output-directory | N/A | Where to write simulation data to |
| calculation-method | "Tree" | Specifies how forces should be calculated |
| accuracy-criterion | "Dynamical" | Specifies how the tree approximation scheme |
| tree-accuracy-paramter | 1e-3 / 0.3 | Specifies how strict the tree's discrimination is |
| timestep-accuracy-parameter | 0.01 | Specifies how strict dynamic timestep is |
| batch-duration | 0.1 | How often data is written to the disk |
| max-time | 30 | When the simulation should stop |

### calculation-method

Can be either "Direct" meaning it will do the pairwise O(N^2) calculation or "Tree" to use the O(N lnN) tree code

### tree-accuracy-criterion & tree-accuracy-parameter

Can either be "Geometric" meaning tree nodes are approximated as single particles if they occuy an angular size on the size less than some angle called `tree-accuracy-parameter` or "Dynamical" meaning tree nodes are approximated as a single particle the inqualiry below holds where $\alpha$ is the `tree-accuracy-parameter`. Default values for `tree-accuracy-parameter` are 1e-3 if tree is dynamical and 0.3 if tree is geometric.

$$\frac{GG * M_{node} * s^2}{d^4} < \alpha a $$

where in this formula s is the size of the node, d is the distance to it and a is the acceleration of the particle being acted on at the previous timestep.

### timestep-accuracy-parameter

The dynamical timestep is determined via the formula below where a,j, and s are accleration, jerk, and snap. Eta is then the `timestep-accuracy-parameter`. This formula is evaluated for each particle and the global timestep is taken to be the minimun value.

$$\Delta t = \frac{\eta}{\sqrt{\left(\frac{j}{a}\right)^2 + \frac{s}{a}}}$$

# [Diagnostics]

The diagnostics keys simply signal which if any simple diagnostic plots to produce. All except the directory paths are simple booleans indicating if the plot should be made or not

| Key  | Description |
|-----|-------------|
| output-directory | Where to write plots to |
| data-directory | Where to find the simulation data |
| trajectories | makes plot of trajectories |
| energy | makes plot checking energy conservation |
| momentum | makes plot checking momentum conservation |
| angular-momentum | makes plot checking angular momentum conservation |