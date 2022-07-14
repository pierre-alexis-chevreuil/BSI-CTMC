function [time_max_Ex, time_max_Ey, maxEx, maxEy, ellipticity_Ex, orientation_Ex, ellipticity_Ey, orientation_Ey] = findPeakElectricField(t, E, time_window)
%% Code to find the parameters of the complex electric field around the
% maximum around the local maxima
%
% Example:
%     angle_HWP           = 10 * pi / 180;
%     angle_QWP           = 79 * pi / 180;
%     lambda_0            = 800e-9;
%     lambda_list         = 500e-9 : 1e-9 : 1000e-9; % A broad range: it should cover the entire spectrum. It will be used to interpolate data
%     tau_FWHM            = 10e-15;
%     thickness_window    = 5e-3;
%     t                   = (-200 : 0.1 : 200) * 1e-15;
%     dispersion          = 1;
%     hwpType             = 'RAC5_2'; % 'ideal', 'RAC5_2' 'RSU1_2', 'RSU2_2' (we have achromatic and superachro range 2 RSU2_2)
%     qwpType             = 'RAC5_4'; % 'ideal', 'RAC5_4' 'RSU1_4', 'RSU2_4' (we have achromatic and superachro range 2 RSU2_4)
%     time_window         = 2.7e-15; % To look at the instantaneous ellipticity for a particular cycle, default: [-2.7e-15   2.7e-15] (1 cycle, 1 cycle after)
%     tau                 = tau_FWHM / (2*sqrt(2*log(2)));
%     IEnv                = exp(-t.^2 ./ (2*tau^2));
%     E                   = zeros(2, length(t));
%     E(1,:)              = sqrt(IEnv) .* exp(1i*(2*pi*299792458/lambda_0 * t));
%     [~, data_UVFS]      = getRefractiveIndex('UVFS', lambda_list);
%     [~, spectral_phase] = getGDCurve(lambda_list, lambda_0, [(-8) * 1e-15, ...
%                                                         (-287)*(1e-15)^2, ...
%                                                         (-236)*(1e-15)^3 ...
%                                                         (-200)*(1e-15)^4]);
%     [t, E]              = chirperPlate(t, E, lambda_list, spectral_phase, dispersion);
%     [t, E]              = quarterWP_PAC(t, E, angle_QWP, qwpType, dispersion);
%     [t, E]              = halfWP_PAC(t, E, angle_HWP, hwpType, dispersion);
%     [t, E]              = chirperPlate(t, E, lambda_list, data_UVFS.spectral_phase_value * thickness_window, dispersion); % Beamline window
%     [t, E]              = removeGDOffset(t, E, lambda_list, lambda_0, dispersion);
%     dispersion_values                           = getDispersionValuesFromEField(t, E(1, :), lambda_0, 3, true);
%     [ellipticity, ExHull, EyHull, orientation]  = ellipticity_finder(t, E, 0, time_window);
%     Ex                                          = real(E(1, :));
%     Ey                                          = real(E(2, :));
%     [time_max_Ex, time_max_Ey, maxEx, maxEy, ellipticity_Ex, orientation_Ex, ellipticity_Ey, orientation_Ey] = findPeakElectricField(t, E, time_window);
%     
%     disp(['Overall ellipticity....... ', num2str(ellipticity)])
%     disp(['Orientation (degrees)..... ', num2str(orientation)])
%     disp(' ')
%     disp('............. E field x .............')
%     disp(['   Max field ............. ', num2str(maxEx)])
%     disp(['   Time max field......... ', num2str(time_max_Ex*1e15), ' fs'])
%     disp(['   Ellipticity max field.. ', num2str(ellipticity_Ex)])
%     disp(['   Orientation max field.. ', num2str(orientation_Ex)])
%     disp(['   Pulse duration......... ', num2str(FWHM(t, abs(E(1, :)).^2) * 1e15), ' fs'])
%     disp(' ')
%     disp('............. E field y .............')
%     disp(['   Max field ............. ', num2str(maxEy)])
%     disp(['   Time max field......... ', num2str(time_max_Ey*1e15), ' fs'])
%     disp(['   Ellipticity max field.. ', num2str(ellipticity_Ey)])
%     disp(['   Orientation max field.. ', num2str(orientation_Ey)])
%     disp(['   Pulse duration......... ', num2str(FWHM(t, abs(E(2, :)).^2) * 1e15), ' fs'])
%     
%     subplot 221
%     plot(t * 1e15, Ex)
%     xlabel('Time (fs)')
%     ylabel('E field (o axis)')
%     grid on
%     ylim([-1 1])
%     subplot 222
%     plot(t * 1e15, Ey)
%     xlabel('Time (fs)')
%     ylabel('E field (e axis)')
%     grid on
%     ylim([-1 1])
%     subplot 223
%     plot(Ex, Ey)
%     hold on
%     plot(ExHull, EyHull, 'b--')
%     hold off
%     xlabel('Time (fs)')
%     xlabel('E field (o axis)')
%     ylabel('E field (e axis)')
%     grid on
%     xlim([-1 1])
%     ylim([-1 1])
%     axis square
%     subplot 224
%     plot3(t * 1e15, Ex, Ey);
%     xlabel('Time (fs)')
%     ylabel('E field (o axis)')
%     zlabel('E field (e axis)')
%     ylim([-1 1])
%     zlim([-1 1])
%     view(3)
%     grid on
%     drawnow
%     makeItNice('FullScreen', false)
%     sgtitle(['Simulation for ', hwpType, ' and ', qwpType, ' waveplates'], 'fontsize', 25, 'interpreter', 'none')
%
% Date: 03.03.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Parameters
Ex                  = real(E(1, :));
Ey                  = real(E(2, :));
[maxEx, indexMaxEx] = max(Ex);
[maxEy, indexMaxEy] = max(Ey);
time_max_Ex         = t(indexMaxEx);
time_max_Ey         = t(indexMaxEy);
if nargout > 4 % Calculates only if needed: computation intensive
    [ellipticity_Ex, ~, ~, orientation_Ex]  = ellipticity_finder(t, E, time_max_Ex, time_window);
    [ellipticity_Ey, ~, ~, orientation_Ey]  = ellipticity_finder(t, E, time_max_Ey, time_window);
end
end