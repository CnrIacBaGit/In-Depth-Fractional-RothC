https://github.com/CnrIacBaGit/In-Depth-Fractional-RothC/genetic

Tested on Ubuntu 22.04, GCC 11.4

Main file:

z_rothC_genetic.py - GA solver for mixed-value optimisation problem for restoration measures selection with process-based simulations performed by z_rothC solver
z_rothC_genetic.json - file with GA solver parameters. Variable with the name in form "__{name}" is a string that describe the meaning of the variable {name}.

Input files (each line - {time, s} {value}):
ET.txt - evapotranspiration, m/s
Kc.txt - soil cover coefficient
SOC_inputs.txt - total SOC inputs, kg m^{-2} s^{-1}
Ta.txt - temperature, oC
gw_depth.txt - water table level, m
precipitation.txt - precipitation, m/s
vgm - layer-based soil model (van Genuchten-Mualem models coefficients, specific storage, thermal coefficients, diffusion coefficients of SOC compounds). Meaning of the values is given in the first line of the file.

Auxiliary scripts:

build.sh - compiles the solver;
sobol_index/sobol_index.py - calculates Sobol indices for 4 variables on the base of GA output (input.csv)
