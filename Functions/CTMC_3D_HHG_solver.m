function [t_CTMC, x_CTMC, vx_CTMC, y_CTMC, vy_CTMC, z_CTMC, vz_CTMC] = CTMC_3D_HHG_solver(time, Ex, Ey, bool_Coulomb_potential, ionization_level, vector_initial_conditions, ODE_rel_tol, ODE_abs_tol, bool_B_field, soft_core_param) %#codegen
%% I use atomic units inside the code (less multiplications needed, and
% easier debugging), but inputs and outputs are in SI, conversion is done
% within the code.
% I tried to optimize the code as much as possible, but it's far from being
% perfect. I use ode78 since it seems to be faster than ode45 for our
% particular set of ODEs. To gain a bit of speed, I also deleted the ion
% potential when the electron is too far away from it.
%
% Important note in May 2022: I added B-field, so electron can be pushed in
% beam propagation direction (along z). I assume the electric field and
% magnetic fields to be 0 along that direction, since the laser should have
% components only along x and y axis
%
% Another important note in June 2022: The potential is now taken WITH the
% soft-core parameter. It is less physical and changes slightly the
% prediction of saddle-point and barrier points, but it lets less chances
% to the electron to start from a point slightly towards the inside of the
% barrier, which can trap it in case of BSI and start at saddle point ...
%
% Inputs:   - time: Time list STARTING AT IONIZATION TIME, not at pulse
%               start
%           - Ex: Real part of E field projected along x axis, at times
%               given in the vector time
%           - Ey: Same as Ex along y
%           - bool_Coulomb_potential: 0 or 1, to account for ion potential
%               or not
%           - ionization_level: average ionization level in the plasm, at
%               times given in the vector time (between 0 and inf, 0 is at
%               pulse beginning. Don't add 1 for ejected electron, it's
%               done in the code!)
%           - vector_initial_conditions: [x_start, vx_start, y_start, vy_start]
%           - ODE_rel_tol: Relative tolerance for the ODE. Can only be one
%               of these: 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8 or 1e-9
%           - ODE_abs_tol: Same as above, but for absolute tolerance
%           - bool_B_field: Include magnetic component of Lorentz force
%               (relevant only for high intensities or long wavelengths),
%               has to be 0 or 1.
%
% Outputs:  - t_CTMC: DIFFERENT FROM TIME VECTOR INPUT !!!!!!!!!!!!!!!!!!!!
%               !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
%           - x_CTMC: distance to core in a.u., at times given in t_CTMC
%           - other outputs: clear
% 
% Date: 22.04.2022 (16.05.2022: added B-field, 13.06.2022: soft-core
% changed)
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

conv_elec_field_au  = 5.1422e11;
conv_time_au        = 2.4189e-17;
conv_au_pos         = 0.5292e-10;
conv_au_speed       = 2.1877e6;
numel_Ex            = numel(Ex);
dt                  = (time(2) - time(1)) / conv_time_au;

if     ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-3 % Matlab doesn't support non-constants inputs for odeset
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-3
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-3
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-3
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-3
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-3
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-3
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-3, 'Stats', 'off');
elseif ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-4
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-4, 'Stats', 'off');
elseif ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-5
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-5, 'Stats', 'off');
elseif ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-6
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-6, 'Stats', 'off');
elseif ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-7
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-7, 'Stats', 'off');
elseif ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-8
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-8, 'Stats', 'off');
elseif ODE_rel_tol == 1e-3 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-9, 'Stats', 'off');
elseif ODE_rel_tol == 1e-4 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-4, 'AbsTol', 1e-9, 'Stats', 'off');
elseif ODE_rel_tol == 1e-5 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-9, 'Stats', 'off');
elseif ODE_rel_tol == 1e-6 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, 'Stats', 'off');
elseif ODE_rel_tol == 1e-7 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-7, 'AbsTol', 1e-9, 'Stats', 'off');
elseif ODE_rel_tol == 1e-8 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-9, 'Stats', 'off');
elseif ODE_rel_tol == 1e-9 && ODE_abs_tol == 1e-9
    opts = odeset('RelTol', 1e-9, 'AbsTol', 1e-9, 'Stats', 'off');
else % In case something is wrong, go back to what gives roughly ok results
    opts = odeset('RelTol', 1e-5, 'AbsTol', 1e-5, 'Stats', 'off');
end

[t_CTMC, xy_result] = ode78(@(t, y) CTMC_fun_2D(t, y, time(1) / conv_time_au, ... % ODE78 approximately 50% faster than ode45 at 3e16 W/cm^ in Kr
                                                      Ex / conv_elec_field_au, ...
                                                      Ey / conv_elec_field_au, ...
                                                      numel_Ex, ...
                                                      ionization_level, ...
                                                      bool_Coulomb_potential, ...
                                                      dt, ...
                                                      bool_B_field, ...
                                                      soft_core_param), ...
                                                      [time(1) time(end)] ./ conv_time_au, ... % If we give discrete values here (we don't let Matlab find the position at which evaluation is done), it's in principle less stable but it is ~10% faster
                                                      [vector_initial_conditions 0 0], ... % 0 0 for z = 0 and v_z = 0 (beam propagation direction)
                                                      opts);
t_CTMC              = t_CTMC * conv_time_au;
x_CTMC              = xy_result(:, 1) * conv_au_pos;
vx_CTMC             = xy_result(:, 2) * conv_au_speed;
y_CTMC              = xy_result(:, 3) * conv_au_pos;
vy_CTMC             = xy_result(:, 4) * conv_au_speed;
z_CTMC              = xy_result(:, 5) * conv_au_pos;
vz_CTMC             = xy_result(:, 6) * conv_au_speed;
    
% subplot 221
%     plot3(t_CTMC / conv_time_au, x_CTMC / conv_au_pos + Ex(1) / conv_elec_field_au * 200, y_CTMC / conv_au_pos + Ey(1) / conv_elec_field_au * 200)
%     hold on
% subplot 222
%     plot(t_CTMC / conv_time_au, x_CTMC / conv_au_pos)
%     hold on
%     grid on
% subplot 223
%     plot(t_CTMC / conv_time_au, y_CTMC / conv_au_pos)
%     hold on
%     grid on
% subplot 224
%     plot(x_CTMC / conv_time_au, y_CTMC / conv_au_pos)
%     hold on
% subplot 221, view(10, 80), grid on, hold off
% subplot 222, hold off
% subplot 223, hold off
% subplot 224, hold off
% makeItNice('fullScreen', false)

function dxdydt = CTMC_fun_2D(t, xy, start_time, Ex, Ey, numel_Ex, ionization_level, bool_Coulomb_potential, dt, bool_B_field, soft_core_param) % For CTMC
    index_current_time  = round((t - start_time) / dt);
    if index_current_time < 1
        index_current_time = 1;
    elseif index_current_time > numel_Ex
        index_current_time = numel_Ex;
    end

    dxdydt      = zeros(6, 1);
    if bool_B_field % Formulas for cross-product of Lorentz force https://en.wikipedia.org/wiki/Lorentz_force (+ and - inverted because of electron negative)
        B_x             = Ey(index_current_time) / 376.730313669; % I assume that B and E are perpendicular (in principle, not true if epsilon_material != 1, but ok approx for gases in infrared), 376 is vacuum impedance
        B_y             = Ex(index_current_time) / 376.730313669; % I assume that B and E are perpendicular (in principle, not true if epsilon_material != 1, but ok approx for gases in infrared), 376 is vacuum impedance
        Lorentz_force_x = - Ex(index_current_time) + xy(6) * B_y; % q(Ex + vy*Bz - vz*By), q = -e = -1 (a.u.)
        Lorentz_force_y = - Ey(index_current_time) - xy(6) * B_x; % q(Ey + vz*Bx - vx*Bz), q = -e = -1 (a.u.)
        Lorentz_force_z = - xy(2) * B_y + xy(4) * B_x; % q(Ez + vx*By - vy*Bx), Ez = 0 & B_z = 0 by assumption
    else
        Lorentz_force_x = - Ex(index_current_time);
        Lorentz_force_y = - Ey(index_current_time);
        Lorentz_force_z = 0; % No electric field along beam propagation direction
    end

    if xy(1)*xy(1) + xy(3)*xy(3) + xy(5)*xy(5) > 400 % if sqrt(x^2 + y^2 + z^2) > 20 a.u. (if outside of Coulomb force)
        dxdydt(1)   = xy(2);
        dxdydt(2)   = Lorentz_force_x;
        dxdydt(3)   = xy(4);
        dxdydt(4)   = Lorentz_force_y;
        dxdydt(5)   = xy(6);
        dxdydt(6)   = Lorentz_force_z;
    else % Around Coulomb
        Z_c         = 1 + ionization_level(index_current_time);
        denominator = soft_core_param + xy(1)*xy(1) + xy(3)*xy(3) + xy(5)*xy(5); % 1 + x^2 + y^2 + z^2
        Coulomb     = - Z_c / realsqrt(denominator * denominator * denominator); % HAS TO BE MINUS, 100% SURE Litteral description since x^2 is slower than x*x
        dxdydt(1)   = xy(2);
        dxdydt(2)   = Lorentz_force_x + bool_Coulomb_potential * Coulomb * xy(1);
        dxdydt(3)   = xy(4);
        dxdydt(4)   = Lorentz_force_y + bool_Coulomb_potential * Coulomb * xy(3);
        dxdydt(5)   = xy(6);
        dxdydt(6)   = Lorentz_force_z + bool_Coulomb_potential * Coulomb * xy(5);
    end
end
end