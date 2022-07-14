function [ellipticity, ExHull, EyHull, orientation_axis_in_degrees] = ellipticity_finder(time, E, center_time, delta_time) %#codegen
%% Code to get the ellipticity of the electric field (typically after some
% waveplates). The ellipticity is given for the electric field, and
% corresponds to the ratio between the real parts of the electric fields
% along the x axis and the y axis. The ellipticity is comprised between -1
% and 1.
%   -1  is left handed circular polarization
%   0   is linear polarization
%   1   is right handed circular polarization
% The ellipticity is evaluated by taking the most outter part of the
% electric field contour (convex Hull). The ellipticity may be different at
% the beginning and the end of the pulse, if the pulses are really short
% because of different GDD along the different axis.
% 
% Input: - time: list of times in seconds
%        - E: complex electric field along x and y (therefore nx2 or 2xn,
%               conversion automatically done in code)
%
% Outputs: - ellipticity: -1 to 1
%          - ExHull: Contour of the convex Hull used for the ellipse fit
%          - EyHull: Same as above, but for y axis
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
% Date: 07.04.2022 (modification of orientation calculation, now based on
%  the ratio of Hull along x and y)
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Parameters
low_thresh              = 0.15; % 0.15 Don't touch: threshold for detecting where the signal on Ex and Ey is non zeros to find orientation of the ellipse (0 to 1)
high_thresh             = 0.85; % 0.85 Don't touch: threshold for detecting where the signal on Ex and Ey is non zeros to find orientation of the ellipse (0 to 1)
orientation_tolerance   = 1e-3;

%% Checking input and resizing data if needed
if size(E, 1) ~= 2 && size(E, 2) ~=2
    error(['The electric field should be of dimension 2xn or nx2 (projections along the x and y axis). You sent a vector of dimension ', num2str(size(E, 1)), ' x ', num2str(size(E, 2)), '. Not good man, not good ...']);
end

if nargin > 1 % Multiple inputs, so user wants to get the ellipticity for a given time range only
    if nargin ~= 4
        error('In this code, you can send either one argument (electric field), either 4 (electric field, time list, center time and delta time), but not less or more. See help of the function to understand why.');
    end
    if numel(delta_time) > 1
        error('Delta time should contain only one element (time considered is from center_time - delta_time to center_time + delta_time)')
    end
    index_start = findIndex(time, center_time - delta_time);
    index_end   = findIndex(time, center_time + delta_time);
    if index_start < 1
        index_start = 1;
    end
    if index_end > numel(time)
        index_end = numel(time);
    end
    if index_start == index_end
        error(['Time entered is incorrect: center time is ', sprintf('%4.2f', center_time*1e15), ' fs and delta time is ', sprintf('%4.2f', delta_time*1e15), ' fs']);
    end
    indices_to_keep = index_start : index_end;
    E               = E(:, indices_to_keep);
end

tempo_amplitude_Ex  = real(abs(E(1, :)));
E(1, :)             = real(E(1, :) ./ tempo_amplitude_Ex); % Renormalization, otherwise ellipticity depends on the difference of electric field between beginning and end of the Hull window
E(2, :)             = real(E(2, :) ./ tempo_amplitude_Ex); % Renormalization, otherwise ellipticity depends on the difference of electric field between beginning and end of the Hull window

%% Getting ellipticity
Ex          = real(E(1, :));
Ey          = real(E(2, :));
if numel(nonzeros(Ex)) < 10 || numel(nonzeros(Ey)) < 10
    ellipticity = 0;
    ExHull                      = Ex;
    EyHull                      = Ey;
    orientation_axis_in_degrees = NaN; % No sense to define 
else
    kHull       = convhull(Ex, Ey);
    ExHull      = Ex(kHull) * max(tempo_amplitude_Ex); % Put back at right amplitude (not important for this code, but it can be externally overlapped with real(E))
    EyHull      = Ey(kHull) * max(tempo_amplitude_Ex); % Put back at right amplitude (not important for this code, but it can be externally overlapped with real(E))
    x           = ExHull' - mean(ExHull, 2);
    y           = EyHull' - mean(EyHull, 2);
    X           = [x.^2, x.*y, y.^2, x, y];
    a           = sum(X) / (X' * X);
    [a,b,c,d,e] = deal(a(1), a(2), a(3), a(4), a(5));
    [~, index_t_max_E_hull]         = max(ExHull.^2 + EyHull.^2);
    orientation_axis_in_degrees     = atand(ExHull(index_t_max_E_hull) / EyHull(index_t_max_E_hull));

    if ( min(abs(b/a), abs(b/c)) > orientation_tolerance)
        orientation_rad = 1/2 * atan( b / (c-a) );
        cos_phi         = cos( orientation_rad );
        sin_phi         = sin( orientation_rad );
        [a,b,c,d,e]     = deal( a*cos_phi^2 - b*cos_phi*sin_phi + c*sin_phi^2,...
                                0,...
                                a*sin_phi^2 + b*cos_phi*sin_phi + c*cos_phi^2,...
                                d*cos_phi - e*sin_phi,...
                                d*sin_phi + e*cos_phi );
    end
    
    result_fit = a*c;

    if (result_fit > 0)
        if (a<0), [a,c,d,e] = deal( -a,-c,-d,-e ); end
        F                   = 1 + (d^2)/(4*a) + (e^2)/(4*c);
        [a,b]               = deal( sqrt( F/a ),sqrt( F/c ) );    
        current_ellipticity = a/b;
    else
        current_ellipticity = 0; % Probably it's simply linear polarization and code doesn't like it
    end

    if current_ellipticity == 1
        orientation_axis_in_degrees = NaN; % meaningless to talk about an orientation for a fully circular polarization
    end
    
    %% Getting direction (left- or right-handed)
    if current_ellipticity > 1, current_ellipticity = 1 / current_ellipticity; end
    indices_possible        = ~((abs(Ex) < low_thresh * max(abs(Ex))) | (abs(Ex) > high_thresh * max(abs(Ex)))) & ...
                              ~((abs(Ey) < low_thresh * max(abs(Ey))) | (abs(Ey) > high_thresh * max(abs(Ey))));
    indices_possible        = nonzeros(double(indices_possible) .* (1:length(indices_possible)));
    index_to_treat          = indices_possible(1);
    if index_to_treat == 1
        index_to_treat = 2;
    end
    x1          = Ex(index_to_treat-1);
    x2          = Ex(index_to_treat);
    x3          = Ex(index_to_treat+1);
    y1          = Ey(index_to_treat-1);
    y2          = Ey(index_to_treat);
    y3          = Ey(index_to_treat+1);
    slope       = (y3-y1) / (x3-x1);
    expected_y2 = slope * x2 + y3 - slope * x3;
    ellipticity = current_ellipticity * round((y2 - expected_y2) ./ abs(expected_y2 - y2)) * round((x1 - x3) ./ abs(x3 - x1));
end
end