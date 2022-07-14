function varargout = getGDCurve(lambda_vector, lambda_central, dispersion_matrix)
%% getGDCurve allows one to get the group delay, the phase (in radians) from
% the dispersion coefficients (GD in s, GDD in s^2, TOD in s^3 ...)
%
% Input arguments : - lambda_vector : Vector containing all the wavelengths
%                       considered (in meters)
%                   - lambda_central : Wavelength around which the
%                       derivatives will be computed (in meters)
%                   - dispersion matrix : Array containing all the
%                       dispersion elements (GD, GDD, TOD, FOD, ... in seconds^n)
% Example :
%
%     lambda      = (655:765) * 1e-9;
%     GD_offset   = -1000     * 1e-15;
%     GDD         = -2323     * 1e-15^2;
%     TOD         = -10365    * 1e-15^3;
%     FOD         = 0         * 1e-15^4;
%     [GD_vs_wavelength , ...
%         phase_vs_wavelength] = getGDCurve(lambda , 730e-9 ,[GD_offset GDD TOD FOD]);
%     yyaxis left
%     plot(lambda * 1e9 , GD_vs_wavelength * 1e12)
%     ylabel('GD (ps)')
%     yyaxis right
%     plot(lambda * 1e9 , phase_vs_wavelength)
%     ylabel('phase (rad)')
%     makeItNice
%
% Author : Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch), adapted from
% Justinas Pupeikis code
%
% Last modification : 05.04.2022 (modified some stuff to make it
% compatible with the compiler)

%% Reminders
% phase(omega) = phase0(omega) + dPhi/dOmega*(omega-omega0) + 1/2!*(omega-omega0)^2*d^2Phi/dOmega^2 + ...
% phase(omega) = phase0(omega) + GD(omega)*(omega-omega0)   + 1/2!*(omega-omega0)^2*GDD + ...
% GD(omega)    = GD            + (omega-omega0) * GDD       + 1/2! * (omega-omega0)^2 * TOD + ...

%% Initialization
c                       = 299792458;
omega_vector            = 2 * pi * c ./ lambda_vector;
omega_central           = 2 * pi * c ./ lambda_central;
GD_vs_lambda            = zeros(1 , length(lambda_vector)); % Initialization of group delay
phase_vs_lambda         = zeros(1 , length(lambda_vector));

%% Computes dispersion
for k = 1 : length(dispersion_matrix) % Reminder for the computation (page 6) : https://ethz.ch/content/dam/ethz/special-interest/phys/quantum-electronics/ultrafast-laser-physics-dam/education/lectures/ultrafast_laser_physics/lecture_notes/3_Dispersion%20compensation.pdf
    phase_vs_lambda = phase_vs_lambda + dispersion_matrix(k) .* (omega_vector - omega_central).^k / factorial(k);
    if abs(dispersion_matrix(k)) > (100 * (1e-12)^k) % Checks if the amount of dispersion added is reasonable or not
        error(['You tried to get add a dispersion of ' , sprintf('%d', int64(dispersion_matrix(k) * 10^(12*k))), ...
               ' fs^' , sprintf('%d', int64(k)), '. That sounds like a lot. Are you sure ? The code expects an input dispersion matrix ', ...
               'under the form [GD GDD TOD FOD ...], eveyrthing in s^n (first is GD in s, not GDD in s^2 !!!).', ...
               'The code will continue normally, but be careful with the result !!']);
    end
end

for k = 1 : length(dispersion_matrix) % Reminder for the computation (page 6) : https://ethz.ch/content/dam/ethz/special-interest/phys/quantum-electronics/ultrafast-laser-physics-dam/education/lectures/ultrafast_laser_physics/lecture_notes/3_Dispersion%20compensation.pdf
    GD_vs_lambda    = GD_vs_lambda + (dispersion_matrix(k).*(omega_vector - omega_central).^(k-1)) / factorial(k-1); % Factorial(0) = factorial(1) = 1
end 

%% Output parser
varargout{1}    = GD_vs_lambda;
varargout{2}    = phase_vs_lambda;
end