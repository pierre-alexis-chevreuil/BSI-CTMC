function [t_out, E_out] = arbitraryWP_PAC(t_in, E_in, alpha, retardance) %#codegen
%% Code to apply a waveplate with an arbitrary retardation to an incoming
% electric field. The code uses Jones formalism for calculations, and uses
% the full phase applied by the waveplate, i.e. it accounts for dispersion.
%
% Input:    - t_in: list of times in seconds
%           - E_in: complex electric field along x and y (therefore nx2 or
%                   2xn (conversion automatically done in code)
%           - alpha: Angle (in radian) of the waveplate
%           - retardance: Retardance in radians (lambda/2 : pi)
%
% Outputs:  - t_out: Time list
%           - E_out: Output complex electric field (2 x n)
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
% Date: 01.04.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch), based on a code
% from Fabian Brunner (F)

c               = 299792458;
rotField        = [cos(alpha), -sin(alpha); sin(alpha), cos(alpha)];
rotFieldBack    = [cos(-alpha), -sin(-alpha); sin(-alpha), cos(-alpha)];
E_rot           = rotField * E_in;
L               = size(E_rot, 2);
timeStep        = t_in(2) - t_in(1);
timePd          = timeStep * (0 : L-1) + t_in(1);
timeSpan        = timePd(end) - timePd(1);
freq            = (-L/2:L/2-1) / timeSpan;
E_freq          = (fftshift(fft(E_rot, [], 2), 2));
phaseShiftWP    = ones(size(freq)) * retardance;
axisWP          = ones(size(freq)) * 0;
E_trans         = complex(zeros(size(E_freq)));

for index = 1:length(freq)
    retM                = [1 0;0 exp(1i * phaseShiftWP(index))];
    angM                = [cos(axisWP(index)) -sin(axisWP(index)); sin(axisWP(index)) cos(axisWP(index))];
    JonesM              = angM' * retM * angM;
    IE                  = (E_freq(:,index));
    out                 = JonesM * IE;
    E_trans(:,index)    = out;
end

L           = size(E_trans, 2);
freqStep    = freq(2) - freq(1);
freqPd      = freqStep*(0:L-1) + freq(1);
freqSpan    = freqPd(end) - freqPd(1);
t_out       = (-L/2:L/2-1) ./ freqSpan; % time vector matching to fftshifted spectrum
E_out       = ifft(ifftshift(E_trans,2), [], 2);

for index = 1 : length(E_out)
    E_out(1:2,index)  = rotFieldBack*[E_out(1,index); E_out(2,index)];
end

checkChirpingFFT(t_out, E_out); % Checks if pulse duration is not too long
end