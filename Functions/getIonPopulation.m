function [list_random_times, ion_origin_electrons_random_times, w_nm, pop_atoms, ionization_level] = ...
    getIonPopulation(t, E, gas, max_ion, PPT_or_ADK, numTrajectories, time_window_ionization, linear_or_circular) %#codegen
%% Code to get the ionization rates, ionic population, and random
% ionization times following the distribution given by the ionization 
% rates.
%
% Inputs:   - t: Time list
%           - E: Electric field (real part only, 1 x n)
%           - gas: 'He', 'Ne', ... as in ADK code
%           - max_ion: 1 for xx+, 2 for xx2+, ...
%           - PPT_or_ADK: 'PPT', 'ADK', 'ADK_BSI', ...
%           - numTrajectories: if some random times of ionization are
%               desired (for CTMC), number of such values. Min value is 1
%           - linear_or_circular: 'lin' or 'circ', for the calculation of
%               the ADK rates
%
% Outputs:  - list_random_times: List of random times (1 x n), with n
%               corresponding to numTrajectories. Those emission times
%               follow the ionization rates distribution
%           - ion_origin_electrons_random_times: At some times, electrons
%               can originate from different ionic states (xx1+, xx2+ ...).
%               This vector gives the origin of the state from which the
%               electron has been emitted (follows the instantaneous
%               ionic population distribution)
%           - w_nm: Ionization rates for the different species
%               (size number_ionic_states x length_time_vector)
%           - pop_atoms: Same as above, but for the ionic population
%               (i.e. accounts for depletion)
%           - ionization_level: Ionization level from 0 to inf. 1 means
%               that there is 0% of xx, 100% of xx+, and 0% of xx2+, or 10%
%               of xx, 80% of xx+, and 10% of xx+. It is the instaneous
%               average ionic charge.
%
% Example:
%
%     tau_FWHM        = 10e-15;
%     t               = linspace(-100e-15, 100e-15, 10000);
%     lambda          = 800e-9;
%     gas             = 'Kr';
%     PPT_or_ADK      = 'ADK_BSI';
%     I_peak          = 2.6e16*1e4;
%     max_ion         = 8;
%     f_t             = functionGenerator(t, tau_FWHM, 0, 'Gaussian'); % To be scaled later in W/m^2
%     E               = sqrt(753.460626923542 * I_peak * f_t).* cos(2*pi*299792458/lambda*t);
%     numTrajectories = 1000;
%     time_window_ion = [-15e-15 15e-15]; % Assumes that ionization rate is 0 (1e-200 not to have issues) outside of this window (makes it faster)
%     linearOrCircular= 'lin';
%     [list_random_times, ion_origin_electrons_random_times, w_nm, pop_atoms, ionization_level] = getIonPopulation(t, E, gas, max_ion, PPT_or_ADK, numTrajectories, time_window_ion, linearOrCircular);
%     leg = {gas};
%     for i = 1 : max_ion
%         leg{end+1} = [gas, num2str(i), '+'];
%     end
%     
%     subplot 221
%         plot(t*1e15, E)
%         xlabel('Time (fs)')
%         ylabel('Electric field (V.m^{-1})')
%         xlim([-tau_FWHM tau_FWHM]*2e15)
%     subplot 222
%         semilogy(t*1e15, w_nm)
%         ylim([1 1e20])
%         legend(leg)
%         xlabel('Time (fs)')
%         ylabel('Ionization rate w (s^{-1})')
%         xlim([-tau_FWHM tau_FWHM]*2e15)
%     subplot 223
%         plot(t*1e15, pop_atoms')
%         ylim([1e-5 1])
%         legend(leg)
%         xlabel('Time (fs)')
%         ylabel('Ion population')
%         xlim([-tau_FWHM tau_FWHM]*2e15)
%     subplot 224
%         plot(t*1e15, ionization_level, 'k', 'HandleVisibility','off')
%         hold on
%         colors = lines(max_ion);
%         for i = 1 : max_ion
%             index_matches_current_ion = ion_origin_electrons_random_times(ion_origin_electrons_random_times == i-1);
%             test = ion_origin_electrons_random_times(ion_origin_electrons_random_times == i-1);
%             plot(list_random_times(ion_origin_electrons_random_times == i-1)*1e15, ion_origin_electrons_random_times(ion_origin_electrons_random_times == i-1) + rand(1, numel(index_matches_current_ion))*0.45, '.', 'color', colors(i, :))
%         end
%         hold off
%         l = legend(leg);
%         title(l, 'Origin from ions')
%         xlabel('Time (fs)')
%         ylabel('Ionization level')
%         xlim([-tau_FWHM tau_FWHM]*2e15)
%
% Date: 17.04.2022 (added switch between linear and circular rates for ADK)
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

    E(E == 0)       = 1e-200; % Precaution
    E(isnan(E))     = 1e-200; % Precaution
    pop_atoms       = ones(max_ion + 1, numel(t)) * 1e-200;
    pop_atoms(1, 1) = 1; % Initial state
    dt              = t(2) - t(1);
    w_nm            = ones(max_ion+ 1, numel(t)) * 1e-200; % Ionization rates in SI
    [~, index_start_ionization] = min(abs(time_window_ionization(1) - t));
    [~, index_end_ionization]   = min(abs(time_window_ionization(2) - t));

    for k = 1 : max_ion + 1
        if k == 1
            name_gas = gas;
        else
            name_gas = [gas, sprintf('%0.0f', k-1)];
        end

        if strcmp(PPT_or_ADK, 'ADK')
%             if ispc
%                 w_nm(k, index_start_ionization : index_end_ionization)  = ADK_fast_mex(E(index_start_ionization : index_end_ionization), name_gas, 0, linear_or_circular);
%             else
                w_nm(k, index_start_ionization : index_end_ionization)  = ADK_fast(E(index_start_ionization : index_end_ionization), name_gas, 0, linear_or_circular);
%             end

        elseif strcmp(PPT_or_ADK, 'ADK_BSI')
%             if ispc
%                 w_nm(k, index_start_ionization : index_end_ionization)  = ADK_fast_mex(E(index_start_ionization : index_end_ionization), name_gas, 1, linear_or_circular);
%             else
                w_nm(k, index_start_ionization : index_end_ionization)  = ADK_fast(E(index_start_ionization : index_end_ionization), name_gas, 1, linear_or_circular);
%             end

        elseif strcmp(PPT_or_ADK, 'PPT')
            error('PPT not yet implemented')

        elseif strcmp(PPT_or_ADK, 'Yudin_Ivanov')
            error('Yudin_Ivanov not yet implemented')

        else
            error('PPT_or ADK can only be "ADK", "ADK_BSI", "PPT" or "Yudin_Ivanov"');

        end
    end

    for current_t = 2 : numel(t)
        pop_atoms(1, current_t)     = pop_atoms(1, current_t-1) * (1 - w_nm(1, current_t-1) * dt); % Ground state doesn't take population anywhere
        for k = 2 : max_ion+ 1
            pop_atoms(k, current_t) = pop_atoms(k, current_t-1) + dt * (w_nm(k-1, current_t-1) * pop_atoms(k-1, current_t-1) - w_nm(k, current_t-1) * pop_atoms(k, current_t-1));
        end
        pop_atoms(pop_atoms(:, current_t) > 1, current_t)    = 1; % Those 3 verifications steps should be done at each time, to avoid propagation of spurious numbers, but not on the full matrix everytime: it's slow!
        pop_atoms(pop_atoms(:, current_t) < 0, current_t)    = 1e-200; % To avoid numerical artefacts
        pop_atoms(isnan(pop_atoms(:, current_t)), current_t) = 1e-200; % To avoid numerical artefacts
    end

    if strcmp(gas, 'He')
        pop_atoms   = [pop_atoms; 1 - sum(pop_atoms, 1)];
        max_ion     = max_ion + 1;
    end
    pop_atoms(1, 1)     = 0;
    ionization_level    = zeros(1, size(pop_atoms, 2));
    for i = 1 : size(ionization_level, 2)
        ionization_level(i) = sum(pop_atoms(:, i) .* (0:max_ion)', 1);
    end

    inst_em_rates                       = gradient(ionization_level);
    inst_em_rates(inst_em_rates < 0)    = 0;
    list_indices_random_times           = randsample(1:numel(ionization_level), numTrajectories, true, inst_em_rates);
    list_random_times                   = t(list_indices_random_times);
    ion_origin_electrons_random_times   = zeros(1, numel(list_random_times));
    grad_population                         = gradient(pop_atoms); % To know when electrons are emitted, we don't wanna have a look at their population, but at their negative population changes (electron emission)!
    grad_population(grad_population >= 0)   = -1e-200; % Only decreasing populations are interesting, because these are the ones to emit electrons for HHG
    rate_emission                           = -grad_population;
%     proba_distribution                      = rate_emission ./ sum(rate_emission, 1);
    for i = 1 : numel(list_random_times) % We know when the electron was emitted, but don't know from which ion it was emitted! We determine it randomly here
        ion_origin_electrons_random_times(i)   = randsample(1:(max_ion+1), 1, true, rate_emission(:, list_indices_random_times(i))) - 1;
    end
end