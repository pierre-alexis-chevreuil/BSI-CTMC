function [t_out, E_out] = chirperPlate(t_in, E_in, lambda_list, spectral_phase, dispersion) %# codegen
%% Code to add some spectral phase on a 2D electric field vector.
%
% Inputs:   - t_in: Time list (1 x n)
%           - E_in: Complex electric field vector (2 x n, 2 being for the
%               2 polarizations).
%           - lambda_list: a list of all lambdas at which spectrum is not
%               zero. It should be rather too broad rather than too narrow.
%               It will be used to interpolate data. For 800, take 500 nm
%               to 1000 nm, it's enough.
%           - spectral_phase: Spectral phase (in radians), with the values
%               corresponding to lambda_list
%           - dispersion: Bool to know if could should be executed or not
%               (useful to engage or disengage dispersion in larger codes)
%
% Outputs:  - t_out: Time list (1 x n)
%           - E_out: Same as E_in, but out, but with some dispersion
%
% Example:
%     angle_HWP           = 10 * pi / 180;
%     angle_QWP           = 79 * pi / 180;
%     angle_parasitic_WP  = 15 * pi / 180;
%     retardance_parasitic_WP = 0 * pi / 180;
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
%     [t, E]              = arbitraryWP_PAC(t, E, angle_parasitic_WP, retardance_parasitic_WP);
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
% Date: 02.03.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

if dispersion
    c           = 299792458;
    L           = size(E_in, 2);
    timeStep    = t_in(2) - t_in(1);
    timePd      = timeStep * (0 : L-1) + t_in(1);
    timeSpan    = timePd(end) - timePd(1);
    freq        = (-L/2:L/2-1) / timeSpan;
    E_freq      = (fftshift(fft(E_in, [], 2), 2));
    phaseShift  = interp1(lambda_list, spectral_phase, c ./ freq, 'linear', 'extrap');
    E_trans     = complex(zeros(size(E_freq)));

    for index = 1 : length(freq)
        JonesM              = exp(-1i*(phaseShift(index))) * [1 0; 0 1];
        IE                  = (E_freq(:, index));
        out                 = JonesM * IE;
        E_trans(:,index)    = out;
    end

    L           = size(E_trans, 2);
    freqStep    = freq(2) - freq(1);
    freqPd      = freqStep*(0:L-1) + freq(1);
    freqSpan    = freqPd(end) - freqPd(1);
    t_out       = (-L/2:L/2-1) ./ freqSpan;
    E_out       = ifft(ifftshift(E_trans, 2), [], 2);

    checkChirpingFFT(t_out, E_out); % Checks if pulse duration is not too long
else
    t_out   = t_in;
    E_out   = E_in;
end
end