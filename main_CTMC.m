% MAIN FILE TO EXECUTE THE CTMC CODE
%
% You can directly execute it, I've put parameters such that one can obtain
% a result from the first execution ("real" results should use different
% setttings of course).
%
% To run this code faster, you can compile several files into binary: one
% can gain one to two orders of magnitude in speed by doing so. The mex
% files can be produced using the coder app from matlab and executing the
% codes in the Examples folder. The files which were designed to be used
% with a compiled version have the mention %#codegen after their
% definition. To switch between mex and non-mex, just change params.use_mex
% from 0 to 1, and the code will execute all the codes with %#codegen with
% their mex version. Important note: the mex files should be named with
% _mex at the end, i.e. if you have code.m, you should have code_mex.mexw64
% for the CTMC code.
%
% Results:  - 1: max distance after 1 cycle (meters)
%           - 2: max distance after 2 cycles (meters)
%           - 3: max distance after 3 cycles (meters)
%           - 4: max energy after 1 cycle (eV), accounts for Ip of ion
%           - 5: max energy after 2 cycles (eV), accounts for Ip of ion
%           - 6: max energy after 3 cycles (eV), accounts for Ip of ion
%           - 7: energy after 1 cycle (eV), accounts for Ip of ion
%           - 8: energy after 2 cycles (eV), accounts for Ip of ion
%           - 9: energy after 3 cycles (eV), accounts for Ip of ion
%           - 10: energy at end of the pulse (eV), accounts for Ip of ion
%           - 11: time of recombination (seconds), taking 0 time at emission
%           - 12: energy of recombination (eV), accounts for Ip of ion
%           - 13: Ip of ion (in a.u.)
%           - 14: BSI or not BSI (bool)
%
% Version 2: Full 3D calculation, accounts for B-field, and external
% control of soft-core Coulomb parameter
% 
% Last update: 13.06.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

addpath(fullfile(pwd, 'Functions'))
ccc
warning('off', 'MATLAB:nearlySingularMatrix') % Just for ellipticity evaluation, no impact on the CTMC
params.angle_HWP                        = 0.05 * pi / 180; % Avoid 0, to prevent having NaNs in the ellipticity calculation
params.angle_QWP                        = 50 * pi / 180; % Avoid 0, to prevent having NaNs in the ellipticity calculation
params.lambda_0                         = 800e-9; % Central wavelength
params.lambda_list                      = 500e-9 : 1e-9 : 1000e-9; % A broad range: it should cover the entire spectrum. It will be used to interpolate data
params.tau_FWHM                         = 10e-15;
params.thickness_window                 = 0e-3; % Beamline window is 5 mm, but best is to keep this 0 and compensate with SLM phase to minimize the number of points in time
params.I_peak                           = 3e16 * 1e4; % W/m^2
params.CEP                              = linspace(0, 2*pi, 1); % CEP randomisation: 50 points is good
params.t                                = (-300 : 0.05 : 300) * 1e-15; % 0.1 fs ok. 0.05 fs very good
params.dispersion                       = 1; % Bool 0 or 1
params.PPT_or_ADK                       = 'ADK_BSI'; % 'ADK', 'ADK_BSI' or 'PPT'
params.hwpType                          = 'RAC5_2'; % 'ideal', 'RAC5_2' 'RSU1_2', 'RSU2_2'
params.qwpType                          = 'RAC5_4'; % 'ideal', 'RAC5_4' 'RSU1_4', 'RSU2_4'
params.time_window                      = 2.7e-15; % To look at the instantaneous ellipticity for a particular cycle, default: 2.7e-15 (1 cycle before and after)
params.recombination_sphere             = 5 * 5.291e-11; % recombination sphere for HHG: 5 atomic radii, following https://iopscience.iop.org/article/10.1088/1361-6455/50/3/035601/ampdf, page 11
params.show_dispersion_values           = false; % true or false, dispersion parameters at HHG spot
params.time_window_ion                  = [-15e-15, 5e-15]; % Assumes that ionization rate is 0 (1e-200 not to have issues) outside of this window (makes it faster), make plot(t, pop_atoms') to see what is realistic (default [-15e-15 5e-15])
params.dispersion_SLM                   = [(-2) * 1e-15, ... % -8 for 5 mm window, -2 for 0 mm
                                           (-106)*(1e-15)^2, ... %-287 for 5 mm window, -106 for 0 mm
                                           (-127)*(1e-15)^3 ... % -236 for 5 mm window, -127 for 0 mm
                                           (-262)*(1e-15)^4]; % -200 for 5 mm window, -262 for 0 mm
params.soft_core_param                  = 0.1; % In a.u. (1 a.u. = 1 Bohr radius, i.e. 0.5e-10 m)
params.time_end_propagation             = 15e-15; % Time after the peak of the pulse after which the code stops propagating the electron
params.grid_size                        = 12; % For getTunnelOrBSIProperties (in atomic units, 8 is default)
params.number_points_grid               = 300; % For getTunnelOrBSIProperties (50 is default for BSI, 100 in tunnel)
params.assumption_tunnel                = 'noAssumption'; % '0speed', '0pos', '0all', 'noAssumption' (don't use 0pos if bool_Coulomb_potential is on: electron will be trapped in potential: long calculation time and useless)
params.bool_Coulomb_potential           = 1; % 0: nothing, 1: soft-core 3D potential with dynamic charge
params.use_mex                          = 0; % 0 or 1. Mex: compiled Matlab code, helps to increase the speed by one to two orders of magnitude. Since it is a binary, a file working on a computer will not necessarily work on another one: you have to compile them yourself one by one
params.include_B_field                  = 1; % 0 or 1
params.plot_tunnel                      = 0; % In the getTunnelOrBSIProperties, should the results be showed (1 only for debug, contains 3D plot, will only work in non-mex mode, don't forget to uncomment last part of the code, commenting is needed for mex compilation)
params.ODE_rel_tol                      = 1e-5; % For the electron trajectory ODE, relative tolerance
params.ODE_abs_tol                      = 1e-5; % For the electron trajectory ODE, absolute tolerance
params.name_mat_results                 = 'Results_CTMC_3D.mat';
params.numTrajectories                  = 1e2; % Per gas and waveplate, with mex and parfor (DON'T FORGET TO PUT BACK PARFOR in get_CTMC_one_HWP_and_QWP_angle if you plot
params.gas                              = 'He';
params.max_ion                          = 1;
params.linear_or_circular               = 'lin';
r                                       = get_CTMC_one_HWP_and_QWP_angle_randomized_CEP(params); % Results

subplot 221
    plot3(r.E_field.t * 1e15, r.E_field.Ex, r.E_field.Ey)
    grid on
    title('Electric field')
    xlabel('Time (fs)')
    ylabel('E-field (v/m)')
    zlabel('E-field (v/m)')
subplot 222
    plot(r.E_field.t * 1e15, r.ion.pop_atoms)
    title('Ionic population')
    xlim([-2 2] * params.tau_FWHM * 1e15)
    xlabel('Time (fs)')
    ylabel('Population')
subplot 223
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 1) * 1e9, 20) % See top of this code to see why it's 1
    xlabel('Distance to origin (nm)')
    ylabel('Number of occurences')
    title('Max distance electron after 1 cycle')
subplot 224
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 10), 20) % See top of this code to see why it's 10
    xlabel('Energy (eV)')
    ylabel('Number of occurences')
    title('Elec kin energy at end of pulse')