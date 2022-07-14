function res_CTMC = get_CTMC_one_HWP_and_QWP_angle(params)
%% Full code to get all parameters of CTMC for one set of HWP and QWP
% waveplates angles (electric field generated, ionization and CTMC).
% This code is made to run fast and calculate tens or hundreds of millions
% of trajectories, so everything has been coded with that in mind. For
% speed reasons, several sub-functions are called through compiled mex64
% files, it therefore requires to run on Windows 64 bits. Don't even try to
% execute this code with more than 10^4-5 points without calling mex, it's
% gonna take forever (speed mex > 100 * matlab code).
% 
% Date: 25.04.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Getting electric field after waveplates
t                   = params.t;
tau                 = params.tau_FWHM / (2*sqrt(2*log(2)));
IEnv                = exp(-t.^2 ./ (2*tau^2));
E                   = zeros(2, length(t));
E(1,:)              = sqrt(753.460626923542 * params.I_peak) * sqrt(IEnv) .* exp(1i*(2 * pi * 299792458 / params.lambda_0 * t + params.CEP));
[~, data_UVFS]      = getRefractiveIndex('UVFS', params.lambda_list, 'needForSpeed', true);
[~, spectral_phase] = getGDCurve(params.lambda_list, params.lambda_0, params.dispersion_SLM);
if params.use_mex
    [t, E]              = chirperPlate_mex(t, E, params.lambda_list, spectral_phase, params.dispersion);
    [t, E]              = removeGDOffset_mex(t, E, params.lambda_list, params.lambda_0, params.dispersion);
    [t, E]              = quarterWP_PAC_mex(t, E, params.angle_QWP, params.qwpType, params.dispersion);
    [t, E]              = halfWP_PAC_mex(t, E, params.angle_HWP, params.hwpType, params.dispersion);
    [t, E]              = chirperPlate_mex(t, E, params.lambda_list, data_UVFS.spectral_phase_value * params.thickness_window, params.dispersion); % Beamline window
    [t, E]              = removeGDOffset_mex(t, E, params.lambda_list, params.lambda_0, params.dispersion);
    [time_max_Ex, time_max_Ey, maxEx, maxEy, ellipticity_Ex, orientation_Ex, ellipticity_Ey, orientation_Ey] = findPeakElectricField(t, E, params.time_window);
else
    [t, E]              = chirperPlate(t, E, params.lambda_list, spectral_phase, params.dispersion);
    [t, E]              = removeGDOffset(t, E, params.lambda_list, params.lambda_0, params.dispersion);
    [t, E]              = quarterWP_PAC(t, E, params.angle_QWP, params.qwpType, params.dispersion);
    [t, E]              = halfWP_PAC(t, E, params.angle_HWP, params.hwpType, params.dispersion);
    [t, E]              = chirperPlate(t, E, params.lambda_list, data_UVFS.spectral_phase_value * params.thickness_window, params.dispersion); % Beamline window
    [t, E]              = removeGDOffset(t, E, params.lambda_list, params.lambda_0, params.dispersion);
    [time_max_Ex, time_max_Ey, maxEx, maxEy, ellipticity_Ex, orientation_Ex, ellipticity_Ey, orientation_Ey] = findPeakElectricField(t, E, params.time_window);
end
Ex                      = real(E(1, :));
Ex(isnan(Ex) | Ex == 0) = 1e-100;
Ey                      = real(E(2, :));
Ey(isnan(Ey) | Ey == 0) = 1e-100;
E_sum                   = sqrt(Ex.^2 + Ey.^2);
res_CTMC.E_field.t                  = t;
res_CTMC.E_field.Ex                 = Ex;
res_CTMC.E_field.Ey                 = Ey;
res_CTMC.E_field.E_sum              = E_sum;
res_CTMC.E_field.max_Ex             = maxEx;
res_CTMC.E_field.max_Ey             = maxEy;
res_CTMC.E_field.max_E_sum          = max(E_sum);
res_CTMC.E_field.time_max_Ex        = time_max_Ex;
res_CTMC.E_field.time_max_Ey        = time_max_Ey;
res_CTMC.E_field.ellipticity_Ex     = ellipticity_Ex;
res_CTMC.E_field.ellipticity_Ey     = ellipticity_Ey;
res_CTMC.E_field.orientation_Ex     = orientation_Ex;
res_CTMC.E_field.orientation_Ey     = orientation_Ey;
res_CTMC.E_field.dispersion_at_HHG  = getDispersionValuesFromEField(t, E(1, :), params.lambda_0, 4, params.show_dispersion_values);

%% Ionization
if params.use_mex
    [list_random_times, ion_origin_electrons_random_times, w_nm, pop_atoms, ionization_level]   = getIonPopulation_mex(t, E_sum, params.gas, params.max_ion, params.PPT_or_ADK, params.numTrajectories, params.time_window_ion, params.linear_or_circular);
    list_Ip_random_times                                                                        = getIonizationPotential(params.gas, ion_origin_electrons_random_times, 'au');  % no mex for this one: nothing computationnaly expensive
    tunnelOrBSIProperties                                                                       = getTunnelOrBSIProperties_mex(t, Ex, Ey, list_random_times, ionization_level, -list_Ip_random_times, params.plot_tunnel, params.grid_size, params.number_points_grid, params.soft_core_param);
else
    [list_random_times, ion_origin_electrons_random_times, w_nm, pop_atoms, ionization_level]   = getIonPopulation(t, E_sum, params.gas, params.max_ion, params.PPT_or_ADK, params.numTrajectories, params.time_window_ion, params.linear_or_circular);
    list_Ip_random_times                                                                        = getIonizationPotential(params.gas, ion_origin_electrons_random_times, 'au');
    tunnelOrBSIProperties                                                                       = getTunnelOrBSIProperties(t, Ex, Ey, list_random_times, ionization_level, -list_Ip_random_times, params.plot_tunnel, params.grid_size, params.number_points_grid, params.soft_core_param);
end
res_CTMC.ion.w_nm                               = w_nm;
res_CTMC.ion.list_random_times                  = list_random_times;
res_CTMC.ion.ion_origin_electrons_random_times  = ion_origin_electrons_random_times;
res_CTMC.ion.pop_atoms                          = pop_atoms;
res_CTMC.ion.ionization_level                   = ionization_level;
res_CTMC.ion.list_Ip_random_times               = list_Ip_random_times;
res_CTMC.ion.tunnelOrBSIProperties              = tunnelOrBSIProperties;

%% CTMC
CTMC_trajectory_analysis  = zeros(params.numTrajectories, 14);
if strcmp(params.assumption_tunnel, '0speed')
    multi_pos_zero_tunnel   = 1;
    multi_speed_zero_tunnel = 0;
elseif strcmp(params.assumption_tunnel, '0pos')
    multi_pos_zero_tunnel   = 0;
    multi_speed_zero_tunnel = 1;
elseif strcmp(params.assumption_tunnel, '0all')
    multi_pos_zero_tunnel   = 0;
    multi_speed_zero_tunnel = 0;
elseif strcmp(params.assumption_tunnel, 'noAssumption')
    multi_pos_zero_tunnel   = 1;
    multi_speed_zero_tunnel = 1;
else % 'noAssumption'
    error('params.assumption_tunnel can only be noAssumption, 0speed, 0pos, or 0all');
end

parfor i = 1 : params.numTrajectories % Solve the equation of motion for each electron and analyzes each trajectory
    [~, index_start_time_list]                  = min(abs(list_random_times(i) - t));
    [~, index_end_time_list]                    = min(abs(params.time_end_propagation - t));
    current_t                                   = t(index_start_time_list : index_end_time_list);
    current_Ex                                  = Ex(index_start_time_list : index_end_time_list);
    current_Ey                                  = Ey(index_start_time_list : index_end_time_list);
    current_list_ion_level                      = ionization_level(index_start_time_list : index_end_time_list);
    if ~isnan(tunnelOrBSIProperties(i, 1))
        if tunnelOrBSIProperties(i, 1) % BSI, initial speed given by difference between potential saddle-point and ionization potential, see getTunnelOrBSIProperties function
            initial_conditions_CTMC = [tunnelOrBSIProperties(i, 2) * multi_pos_zero_tunnel, ...
                                       tunnelOrBSIProperties(i, 8) * multi_speed_zero_tunnel, ...
                                       tunnelOrBSIProperties(i, 3) * multi_pos_zero_tunnel, ...
                                       tunnelOrBSIProperties(i, 9) * multi_speed_zero_tunnel]; % v, vx, y, vy
        else % tunnel
            initial_conditions_CTMC = [tunnelOrBSIProperties(i, 14) * multi_pos_zero_tunnel, ...
                                       0, ... % if tunnel in getTunnelOrBSIProperties, it's a NaN
                                       tunnelOrBSIProperties(i, 15) * multi_pos_zero_tunnel, ...
                                       0]; % v, vx, y, vy
        end

        if params.use_mex
            [t_CTMC, x_CTMC_curr, vx_CTMC_curr, y_CTMC_curr, vy_CTMC_curr, z_CTMC_curr, vz_CTMC_curr]   = CTMC_3D_HHG_solver_mex(current_t, current_Ex, current_Ey, params.bool_Coulomb_potential, current_list_ion_level, initial_conditions_CTMC, params.ODE_rel_tol, params.ODE_abs_tol, params.include_B_field, params.soft_core_param);
            CTMC_trajectory_analysis(i, :)                                                              = findHHGTimesCTMC_mex(t_CTMC, x_CTMC_curr, vx_CTMC_curr, y_CTMC_curr, vy_CTMC_curr, z_CTMC_curr, vz_CTMC_curr, params.recombination_sphere, list_Ip_random_times(i), params.lambda_0, tunnelOrBSIProperties(i, 1));
        else
            [t_CTMC, x_CTMC_curr, vx_CTMC_curr, y_CTMC_curr, vy_CTMC_curr, z_CTMC_curr, vz_CTMC_curr]   = CTMC_3D_HHG_solver(current_t, current_Ex, current_Ey, params.bool_Coulomb_potential, current_list_ion_level, initial_conditions_CTMC, params.ODE_rel_tol, params.ODE_abs_tol, params.include_B_field, params.soft_core_param);
            CTMC_trajectory_analysis(i, :)                                                              = findHHGTimesCTMC(t_CTMC, x_CTMC_curr, vx_CTMC_curr, y_CTMC_curr, vy_CTMC_curr, z_CTMC_curr, vz_CTMC_curr, params.recombination_sphere, list_Ip_random_times(i), params.lambda_0, tunnelOrBSIProperties(i, 1));
        end
%         t_CTMC = t_CTMC - t_CTMC(1); % To plot individual electron trajectories (and change parfor by for)
%         plot(t_CTMC, sqrt(x_CTMC_curr.^2 + y_CTMC_curr.^2 + z_CTMC_curr.^2) / 0.5292e-10, 'color', [0 0.447 0.741])
%         hold on
%         grid on
%         xlim([-0.1e-15 2.7e-15])
%         ylim([0 10])
%         drawnow
    end
end
CTMC_trajectory_analysis(CTMC_trajectory_analysis(:, 8) == 0, :)    = []; % Removes meaningless points
res_CTMC.CTMC.CTMC_trajectory_analysis                              = CTMC_trajectory_analysis;
end