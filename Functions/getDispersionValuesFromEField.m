function dispersion_values = getDispersionValuesFromEField(time, E_field, lambda_0, order_fit, display_values)
%% Code to find the dispersion parameters (GD, GDD, TOD ...) from a
% time-dependent electric field.
%
% Inputs:   - time: time list, one dimension
%           - E_field: Complex electric field versus time (1 dimension: 
%               for different polarizations, you have to project)
%
% Outputs: dispersion_values: [GD (s), GDD (s^2), ...]
%
% Example:
%     angle_HWP           = 80 * pi / 180;
%     angle_QWP           = 120 * pi / 180;
%     lambda_0            = 800e-9;
%     lambda_list         = 500e-9 : 1e-9 : 1000e-9; % A broad range: it should cover the entire spectrum. It will be used to interpolate data
%     tau_FWHM            = 10e-15;
%     thickness_window    = 5e-3;
%     t                   = (-200 : 0.1 : 200) * 1e-15;
%     dispersion          = 1;
%     hwpType             = 'RAC5_2'; % 'ideal', 'RAC5_2' 'RSU1_2', 'RSU2_2' (we have achromatic and superachro range 2 RSU2_2)
%     qwpType             = 'RAC5_4'; % 'ideal', 'RAC5_4' 'RSU1_4', 'RSU2_4' (we have achromatic and superachro range 2 RSU2_4)
%     time_window         = 2 * 2.7e-15; % To look at the instantaneous ellipticity for a particular cycle, default: [-2.7e-15   2.7e-15] (1 cycle, 1 cycle after)
%     tau                 = tau_FWHM / (2*sqrt(2*log(2)));
%     IEnv                = exp(-t.^2 ./ (2*tau^2));
%     E                   = zeros(2, length(t));
%     E(1,:)              = sqrt(IEnv) .* exp(1i*(2*pi*299792458/lambda_0 * t));
%     [~, data_UVFS]      = getRefractiveIndex('UVFS', lambda_list);
%     [~, spectral_phase] = getGDCurve(lambda_list, lambda_0, [(-8) * 1e-15, ...
%                                                             (-287)*(1e-15)^2, ...
%                                                             (-236)*(1e-15)^3 ...
%                                                             (-200)*(1e-15)^4]);
%     [t, E]              = chirperPlate(t, E, lambda_list, spectral_phase, dispersion);
%     [t, E]              = quarterWP_PAC(t, E, angle_QWP, qwpType, dispersion);
%     [t, E]              = halfWP_PAC(t, E, angle_HWP, hwpType, dispersion);
%     [t, E]              = chirperPlate(t, E, lambda_list, data_UVFS.spectral_phase_value * thickness_window, dispersion); % Beamline window
%     [t, E]              = removeGDOffset(t, E, lambda_list, lambda_0, dispersion);
%     dispersion_values               = getDispersionValuesFromEField(t, E(1, :), lambda_0, 3, true);
%     [ellipticity, ExHull, EyHull]   = ellipticity_finder(t, E, 0, time_window);
%     Ex                              = real(E(1, :));
%     Ey                              = real(E(2, :));
%     disp(['Ellipticity............. ', num2str(ellipticity)])
%     disp(['Pulse duration x axis... ', num2str(FWHM(t, abs(E(1, :)).^2) * 1e15), ' fs'])
%     disp(['Pulse duration y axis... ', num2str(FWHM(t, abs(E(2, :)).^2) * 1e15), ' fs'])
%     
%     subplot 221
%         plot(t * 1e15, Ex)
%         xlabel('Time (fs)')
%         ylabel('E field (o axis)')
%         grid on
%         ylim([-1 1])
%     subplot 222
%         plot(t * 1e15, Ey)
%         xlabel('Time (fs)')
%         ylabel('E field (e axis)')
%         grid on
%         ylim([-1 1])
%     subplot 223
%         plot(Ex, Ey)
%         hold on
%         plot(ExHull, EyHull, 'b--')
%         hold off
%         xlabel('Time (fs)')
%         xlabel('E field (o axis)')
%         ylabel('E field (o axis)')
%         grid on
%         xlim([-1 1])
%         ylim([-1 1])
%         axis square
%     subplot 224
%         plot3(t * 1e15, Ex, Ey);
%         xlabel('Time (fs)')
%         ylabel('E field (o axis)')
%         zlabel('E field (e axis)')
%         ylim([-1 1])
%         zlim([-1 1])
%         view(3)
%         grid on
%         drawnow
%     makeItNice('FullScreen', false)
%     sgtitle(['Simulation for ', hwpType, ' and ', qwpType, ' waveplates'], 'fontsize', 25, 'interpreter', 'none')
%
% Date: 02.03.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Inputs checks
if size(time, 1) > size(time, 2)
    time = time';
end
if size(E_field, 1) > size(E_field, 2)
    E_field = E_field';
end
if size(E_field, 1) > 1
    error('This code only accepts 1 x n electric fields. If you sent 2 x n (two polarizations), remove one of the 2.')
end

%% Code
c = 299792458;
L                   = numel(E_field);
timeStep            = time(2) - time(1);
timePd              = timeStep * (0 : L-1) + time(1);
timeSpan            = timePd(end) - timePd(1);
freq                = (-L/2:L/2-1) / timeSpan;
E_freq              = fftshift(fft(fftshift(E_field)));
df                  = 1 / FWHM(time, real(E_field).^2);
f_0                 = c / lambda_0;
spectral_phase      = unwrap(angle(E_freq));
[~, index_min_freq] = min(abs(f_0 - df - freq));
[~, index_max_freq] = min(abs(f_0 + df - freq));
dispersion_values   = getDispersionValues(spectral_phase(index_min_freq: index_max_freq) , ...
                                            c ./ freq(index_min_freq: index_max_freq) , 'phase' , lambda_0, order_fit, ...
                                            'displayValues', display_values);

end