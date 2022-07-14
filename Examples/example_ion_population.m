addpath(fullfile(fileparts(pwd), 'Functions'))
ccc

warning('off', 'MATLAB:legend:IgnoringExtraEntries') % Not happy when plotting origin of ions: nothing emitted from the last ion, since it's not being depleted

tau_FWHM        = 10e-15;
t               = linspace(-100e-15, 100e-15, 10000);
lambda          = 800e-9;
gas             = 'Ar';
PPT_or_ADK      = 'ADK_BSI';
I_peak          = 2.6e16*1e4;
max_ion         = 7;
f_t             = functionGenerator(t, tau_FWHM, 0, 'Gaussian'); % To be scaled later in W/m^2
E               = sqrt(753.460626923542 * I_peak * f_t).* cos(2*pi*299792458/lambda*t);
numTrajectories = 1e2;
time_window_ion = [-15e-15 15e-15]; % Assumes that ionization rate is 0 (1e-200 not to have issues) outside of this window (makes it faster)
xlims           = [-15 10];
linearOrCircular= 'lin';
type_plot       = 1;

[list_random_times, ion_origin_electrons_random_times, w_nm, pop_atoms, ionization_level]   = getIonPopulation(t, E, gas, max_ion, PPT_or_ADK, numTrajectories, time_window_ion, linearOrCircular);

mat_BSI_vs_ions                                                                             = getBSIIntensities(t, I_peak*f_t, gas, max_ion);
leg = {gas};
for i = 1 : max_ion
    leg{end+1} = [gas, num2str(i), '+'];
end

subplot 221
    plot(t*1e15, E)
    xlabel('Time (fs)')
    ylabel('Electric field (V.m^{-1})')
    xlim(xlims)
subplot 222
    semilogy(t*1e15, w_nm)
    ylim([1 1e20])
    legend(leg)
    xlabel('Time (fs)')
    ylabel('Ionization rate w (s^{-1})')
    xlim(xlims)
subplot 223
    colors = lines(max_ion + 1);
    for i = 1 : max_ion + 1
        plot(t*1e15, pop_atoms(i, :), 'color', colors(i, :))
        hold on
        plot(mat_BSI_vs_ions(i, 3) * 1e15 * [1 1], [0 1], '--', 'color', colors(i, :), 'HandleVisibility','off')
    end
    hold off
    ylim([1e-5 1])
    legend(leg)
    xlabel('Time (fs)')
    ylabel('Ion population')
    xlim(xlims)
subplot 224
    plot(t*1e15, ionization_level, 'k', 'HandleVisibility','off')
    hold on
    colors = lines(max_ion);
    for i = 1 : max_ion
        index_matches_current_ion = ion_origin_electrons_random_times(ion_origin_electrons_random_times == i-1);
%         test = ion_origin_electrons_random_times(ion_origin_electrons_random_times == i-1);
        plot(list_random_times(ion_origin_electrons_random_times == i-1)*1e15, ion_origin_electrons_random_times(ion_origin_electrons_random_times == i-1) + rand(1, numel(index_matches_current_ion))*0.45, '.', 'color', colors(i, :))
    end
    hold off
    l = legend(leg);
    title(l, 'Origin from ions')
    xlabel('Time (fs)')
    ylabel('Ionization level')
    xlim(xlims)