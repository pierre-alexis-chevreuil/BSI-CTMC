addpath(fullfile(fileparts(pwd), 'Functions'))
ccc

warning('off', 'MATLAB:nearlySingularMatrix')

%% Parameters
params.angle_HWP                        = 0 * pi / 180;
params.angle_QWP                        = 50 * pi / 180;
params.lambda_0                         = 800e-9;
params.lambda_list                      = 500e-9 : 1e-9 : 1000e-9; % A broad range: it should cover the entire spectrum. It will be used to interpolate data
params.tau_FWHM                         = 10e-15;
params.thickness_window                 = 0e-3; % Beamline window is 5 mm, but best is to keep this 0 and compensate with SLM phase to minimize the number of points in time
params.I_peak                           = 3e16 * 1e4; % W/m^2 (Cutoff Kr @ 1e14 W/cm^2 at 800 nm: 32.9 eV, Ip is 14 eV (32.9-14 = 18.9 = Up max))
params.CEP                              = 0;
params.t                                = (-300 : 0.1 : 300) * 1e-15; % 0.1 fs ok ONLY IF INTERPOLATION FOR ELECTRON DETECTION DONE AFTERWARDS (some electrons can traverse 5 a.u. radius circle within 0.01 fs, i.e. 10 as)
params.dispersion                       = 1; % Bool 0 or 1
params.gas                              = 'He';
params.max_ion                          = 1;
params.numTrajectories                  = 10; % in BSI, ~ 4000 points per sec for getTunnelOrBSIProperties
params.PPT_or_ADK                       = 'ADK_BSI'; % 'ADK', 'ADK_BSI' or 'PPT'
params.hwpType                          = 'RAC5_2'; % 'ideal', 'RAC5_2' 'RSU1_2', 'RSU2_2' (we have achromatic and superachro range 2 RSU2_2)
params.qwpType                          = 'RAC5_4'; % 'ideal', 'RAC5_4' 'RSU1_4', 'RSU2_4' (we have achromatic and superachro range 2 RSU2_4)
params.time_window                      = 2.7e-15; % To look at the instantaneous ellipticity for a particular cycle, default: 2.7e-15 (1 cycle before and after)
params.recombination_sphere             = 5 * 5.291e-11; % 5 atomic radii, following https://iopscience.iop.org/article/10.1088/1361-6455/50/3/035601/ampdf, page 11
params.show_dispersion_values           = false;
params.time_window_ion                  = [-15e-15, 5e-15]; % Assumes that ionization rate is 0 (1e-200 not to have issues) outside of this window (makes it faster), make plot(t, pop_atoms') to see what is realistic (default [-15e-15 5e-15])
params.dispersion_SLM                   = [(-2) * 1e-15, ... % -8 for 5 mm window, -2 for 0 mm
                                           (-106)*(1e-15)^2, ... %-287 for 5 mm window, -111 for 0 mm
                                           (-127)*(1e-15)^3 ... % -236 for 5 mm window, -62 for 0 mm
                                           (-262)*(1e-15)^4]; % -200 for 5 mm window, 250 for 0 mm
params.soft_core_param                  = 1; % In a.u. (1 a.u. = 1 Bohr radius, i.e. 0.5e-10 m)
params.time_end_propagation             = 15e-15; % Time after the peak of the pulse after which the code stops propagating the electron
params.grid_size                        = 8; % For getTunnelOrBSIProperties (in atomic units, 8 is default)
params.number_points_grid               = 150; % For getTunnelOrBSIProperties (50 is default for BSI, 100 in tunnel)
params.assumption_tunnel                = 'noAssumption'; % '0speed', '0pos', '0all', 'noAssumption' (don't use 0pos if bool_Coulomb_potential is on: electron will be trapped in potential: long calculation time and useless)
params.linear_or_circular               = 'lin';
params.bool_Coulomb_potential           = 1; % 0: nothing, 1: soft-core 3D potential with dynamic charge
params.use_mex                          = 0;
params.include_B_field                  = 1;
params.plot_tunnel                      = 0; % In the getTunnelOrBSIProperties, should the results be showed (1 only for debug, contains 3D plot, will only work in non-mex mode)
params.ODE_rel_tol                      = 1e-5; % For the electron trajectory ODE, relative tolerance
params.ODE_abs_tol                      = 1e-5; % For the electron trajectory ODE, absolute tolerance
time_lims                               = [-15 15]; % Just for display in this code

%% Code
tic
r = get_CTMC_one_HWP_and_QWP_angle(params);
toc

%% Plotting
pourcentage_recombination = numel(r.CTMC.CTMC_trajectory_analysis( ~isnan(r.CTMC.CTMC_trajectory_analysis(:, 9)), 9)) ./ numel(r.CTMC.CTMC_trajectory_analysis(:, 9));
subplot 331
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 1) * 1e9, 100)
    title(['Max distance over the first cycle, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 1)*1e9), '%4.2f'), 'nm'])
subplot 332
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 2) * 1e9, 100)
    title(['Max distance over the second cycle, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 2)*1e9), '%4.2f'), 'nm'])
subplot 333
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 3) * 1e9, 100)
    title(['Max distance over the third cycle, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 3)*1e9), '%4.2f'), 'nm'])
subplot 334
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 4), 100)
    title(['Max energy over the first cycle, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 4)), '%4.2f'), 'eV'])
subplot 335
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 5), 100)
    title(['Max energy over the second cycle, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 5)), '%4.2f'), 'eV'])
subplot 336
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 6), 100)
    title(['Max energy over the third cycle, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 6)), '%4.2f'), 'eV'])
subplot 337
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 7), 100)
    title(['Energy at end of the pulse, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis(:, 7)), '%4.2f'), 'eV'])
subplot 338
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 8), 100)
    title(['Time of recombination, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis( ~isnan(r.CTMC.CTMC_trajectory_analysis(:, 8)), 8))*1e15, '%4.2f'), 'fs (', num2str(pourcentage_recombination*100, '%4.2f'), '% recombine)'])
subplot 339
    histogram(r.CTMC.CTMC_trajectory_analysis(:, 9), 100)
    title(['Energy photons after recombination, mean: ', num2str(mean(r.CTMC.CTMC_trajectory_analysis( ~isnan(r.CTMC.CTMC_trajectory_analysis(:, 9)), 9)), '%4.2f'), 'eV'])