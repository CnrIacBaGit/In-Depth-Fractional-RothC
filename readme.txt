Tested on Ubuntu 22.04, GCC 11.4

Main file:

z_rothC [<parameter name> <parameter value>]... - solver for in-depth fractional-order RothC-based SOC dynamics model
    Mandatory parameters - names of input files:
	"vgm" - the name of a file with the van Genuchten-Muallem model's coefficients and other hydrophysical parameters of the soil;
	"et_file" - the name of a file with evapotranspiration values, m/s;
	"prec_file" - the name of a file with precipitation values, m/s;
	"gw_file" - the name of a file with groundwater levels, m;
	"Ta_file" - the name of a file with air temperature values, °C;
	"SOC_file" - the name of a file with the values of SOC inputs, kg/m^2/s;
	"Kc_file" - the name of a file with soil cover coefficient values (k_c in RothC model);
	Each row of "{et,prec,gw,Ta,SOC,Kc}_file" must contain a pair of (<time, seconds>,<value>). Files must be sorted upon time. In simulations values are linearly interpolated.
	In "vgm" file the first row is ignored and contains the names of parameters, other rows represent hydrophysical parameters of soil layers and have the following form:
	    theta_r - residual moisture content (the van Genuchten model);
	    theta_s - saturated moisture content (the van Genuchten model);
	    n – the van Genuchten model's parameter;
	    a – the van Genuchten model's parameter, 100m^{-1};
	    h0,h1 - upper and lower depth of soil layers, m. File must be sorted upon the increase of depth;
	    kf - filtration coefficient, m/s;
	    beta - exponent in the Muallem model;
	    specific_storage - specific storage of soil, m^{-1};
	    lamdba - thermal conductivity of soil, W m^{-1} K^{-1};
	    Ct - heat capacity of soil, J m^{-3} K^{-1};
	    D[0,1,2,3] - diffusivities of SOC compounds in soil, m^2/s;
    Model parameters:
	"C_q" - time-fractional derivative's order q;
	"corrector_type" - the form of dimensions correction factor for fractional derivative (0 - equal to 1, 1 - equal to Gamma(2-q)t^{q-1});
	"fanox" - respiration rate factor for saturated soil (f_anox);
	"inputs_exponent" - exponent of SOC inputs decay with depth;
	"inputs_divisor" - a part of SOC inputs distributed within soil column (1-input_divisor of SOC inputs enters onto surface);
	"oxidation" - a level of CH4 oxidation;
	"ox_tau" - exponent of CH4 oxidation level decay with depth;
	"gamma" - the coefficient of SOC inputs allocation between DPM and RPM pools;
	"C0{0,1,2,3}" - initial value of concentrations for 4 SOC pools;
	"H0" - H0+z is the initial value of water head distribution in iteration process used to determine steady-state solution of Richards equation;
	"soc_inputs_divisor" - a constant, on which SOC inputs are divided during the simulation;
    Simulation parameters:
	"Tm" - ending time in days;
	"NB" - number of cells in a finite-difference grid;
	"LB" - domain depth in meters;
	"testing" - mode 1 - solves artificially constructed testing problem, mode 2 - solves advection-diffusion-reaction equations for SOC components without advection term;
	"fr_eps" - accuracy threshold for fixed memory computation of fractional derivatives;
	"summing_depth" - depth, m, down to which output variables are integrated;
	"boundary_v" - first-order (1) or second-order (2) approximation of flow velocity and its derivative in boundary points;
    Output parameters:
	"Om" - output time interval in days;
	"Zoutstep" - values in cells i*Zoutstep, i=0,...,  will be written in output files;
	"debug_level" - level 1 - outputs info about Pickard iteration process, level 2 - outputs info about linear systems solution;
    Parameters of TFQMR iterative linear systems solver:
	"ls_eps" - accuracy threshold (average sum of squares of differences);
	"ls_max_iter" - maximal number of iterations;
    Parameters of the procedure for dynamic time step change:
	"ls_min_tau" - minimal time step length in seconds (solution fails when this value is reached);
	"max_tau" - upper limit for time step length in seconds;
	"ls_mult" – factor, on which time step length is multiplied or divided when time step changes. Time step length is divided by ls_mult when ls_max_iter iterations are reached during linear system solution and solution procedure is repeated;
	"ls_percent" - if number of TFQMR iteration is lower than ls_percent*ls_max_iter, time step length is multiplied of ls_mult; 
    Output files, in which each row corresponds to a specific moment of time:
	out_H.txt - results of Richards equation solution:
	    t(days) <moment of time in days> tau(seconds) <the value of time step length in seconds> ET <evapotranspiration, m/s> Prec <precipitation, m/s> Fl <the value of water head function on the bottom of the domain, m> - H: <list of water head values, m, in grid cells (values of depth z in the first row)> water_table <simulated water table level, m>; 
	out_T.txt - results of heat transport equation solution:
	    t(days) <moment of time in days> tau(seconds) <the value of time step length in seconds> T: <list of soil temperature values, °C, in grid cells (values of depth z in the first row)>;
	out_C[0,1,2,3].txt - results for SOC compounds dynamic:
	    t(days) <moment of time in days> tau(seconds) <the value of time step length in seconds> SoC input <the total SOC inputs amount, kg/m^2/s> C[{0,1,2,3}]: <list of SOC compound concentration values, kg/m^3, in grid cells (values of depth z in the first row)>;
	    The file out_C3.txt contains in each row additional columns with the values of the following aggregate indicators:
	     CO2 respiration - <moment CO2 respiration from soil column down to summing depth, kg/m^2/s^q> CH4 respiration - <moment CH4 respiration from soil column down to summing depth, kg/m^2/s^q> Balance <moment SOC balance in soil column down to summing depth, kg/m^2/s^q> sumSOC <total (down to summing depth) simulated SOC content in soil column down to summing depth, kg/m^2> sumSoCbalance <total (down to summing depth) SOC content according to balance equation, kg/m^2> sumI <total (down to summing depth) SOC inputs up to the current moment of time, kg/m^2> sumRC02 <total (down to summing depth) CO2 respiration up to the current moment of time, kg/m^2> sumRCH4 <total (down to summing depth) CH4 respiration up to the current moment of time, kg/m^2> sumOut <total SOC compounds outflow through the bottom boundary (summing depth) up to the current moment of time, kg/m^2>;
	out_{V,co2,ch4,rho,soc}.txt - water movement velocity, m/s (V); CO2 respiration, kg/m^3/s^q (co2); CH4 respiration, kg/m^3/s^q (ch4); respiration rate function (rho); total SOC content, kg/m^3 (soc):
	    t(days) <moment of time in days> {V,co2_respiration,ch4_respiration,rho,soc}: <list of values in grid cells (values of depth z in the first row)>;

Auxiliary scripts:

build.sh - compiles the solver;
run_calcs.sh - runs several simulations scenarios, input data for which is stored in ET2.txt, Kc1.txt, Ta1.txt, gw_depth1.txt, vgm{1,2,3}.txt, SOC_inputs2.txt, precipitation2.txt;
move_to_folder.sh <folder name> - moves output files into the specified folder and postprocess them;
transpose.sh <file> - transposes the specified file considered as a matrix with space as a delimiter, writes results to stdout;
yearly_avg.sh <file> - computes yearly sum CO2 and CH4 respiration from the specified file of out_C3.txt format, writes results to stdout; 

