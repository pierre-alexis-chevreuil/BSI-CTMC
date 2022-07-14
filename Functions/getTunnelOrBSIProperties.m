function mat_results = getTunnelOrBSIProperties(t, Ex, Ey, instants_ionization, ionization_level_vs_time, list_Ip, plot_result, grid_size, number_points_grid, soft_core_param) %#codegen
%% Code to get the properties of tunnel or BSI ionization (ionization
% potential, points of exits (x and y) in case of tunneling, potential
% saddle point in case of BSI (along x and y), velocities of exit (0 in
% case of tunneling, not zero for BSI), determines if for the current ion
% we're below or above BSI.
% Note that in case of tunneling, the position of barrier entrance is not
% very precise because the potential gradient is quite steep there.
% However, since we're far from the Coulomb potential at the barrier exit,
% the potential gradient is lower, leading to a descent evaluation of the
% barrier exit position, even without too many points. It's actually quite
% convinient, since the barrier exit position is the most interesting
% parameter (at least for CTMC, since it's the position from which
% propagation is started).
% Since a wide range of intensities can be used here, and that both
% tunneling and BSI should be ok, the grid size is not fixed. What the code
% does is that is starts with a "standard" grid size (typically 16x16 a.u.)
% and will increase either the grid size and keeping the same number of
% points in case of tunneling where tunnel exit is not found, or increase
% the resolution and keeping grid size the same in case the resolution
% became so poor that the soft-core Coulomb potential goes above the
% barrier (in case of tunneling, this is meaningless).
% Typical execution times in mex version (obvisouly without plot):
%   ~ 1000 points per second in tunnel mode
%   ~ 4000 points per second in BSI mode
%
% Inputs:   - t: Vector of time points (only used to find index of
%                   ionization in electric field vectors)
%           - Ex: Electric field vector (real part of E-field along x), in
%                   SI (V/m)
%           - Ey: Same as above for y-axis
%           - instants_ionization: Instants of ionization for each electron
%                   to be studied. It should be one of the elements of t.
%           - ionization_level_vs_time: average ionization level for each
%                   point in time (between 0 and inf)
%           - list_Ip: ionization potential of each atom/ion to be studied
%                   (in atomic units). For exemple, if we have electrons
%                   coming from Kr, Kr2 and Kr4, it's going to be
%                   [0.514, 1.3170, 2.37731]
%           - plot_result: bool (0 or 1) to plot 3D visualization of the
%                   potential laser + ion
%           - grid_size: Size of the 2D-grid to be used for the potential.
%                   It's in atomic units, so 1 a.u. = 5.3e-11 meters.
%                   Default is 16 a.u. (if you want 64 a.u. x 64 a.u. for
%                   example, just put 64, the grid has to be square). If
%                   you don't know what to use, just put 0 and the code
%                   will adapt itself.
%           - number_points_grid: Number of points for the 2D-grid. As for
%                   grid-size, if you don't know what to use, just put 0
%                   and the code will adapt itself.
%
% Output: One struct, with each index being an electron and each name being
% a parameter:
%           1: BSI (1 if BSI, 0 if tunnel)
%           2: x_saddle (x position saddle point, i.e. top barrier)
%           3: y_saddle (y position saddle point, i.e. top barrier)
%           4: r_saddle (sqrt(x_saddle^2 + y_saddle^2))
%           5: height_saddle (height saddle point, i.e. top barrier, not
%                   considering Ip)
%           6: diff_saddle_Ip (difference height saddle point - ionization
%                   potential. If negative, we're in BSI, otherwise we're
%                   in tunnel)
%           7: velocity_after_ionization (in tunnel mode, equals 0. In
%                   BSI mode, corresponds to the diff_saddle_Ip energy
%                   transferred in velocity)
%           8: velocity_after_ionization_x (same as above, but projected
%                   along the x axis, the projection ratio corresponding
%                   to the field strength since force proportional to E)
%           9: velocity_after_ionization_y (same as above for y)
%           10: radiusEntryTunnel (radius of entry of tunnel, in case of
%                   tunneling. In BSI, it's NaN)
%           11: radiusExitTunnel (same as entry, but exit)
%           12: xEntryTunnel (same as above, but natural coordinates)
%           13: yEntryTunnel (same as above)
%           14: xExitTunnel (same as above, but exit)
%           15: yExitTunnel (same as above)
%           16: potentialEntryTunnel (laser + ion potential, no Ip, at
%                   barrier entrance. In BSI, it's NaN)
%           17: potentialExitTunnel (same as above, but exit)
%
% This code is meant to be used with the rest of the CTMC series of codes.
%
% Date: 24.04.2022 (13.06.2022: update soft-core)
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Parameters
if grid_size == 0 && number_points_grid == 0
    if max(sqrt(Ex.^2 + Ey.^2)) < 5e10 % V/m, i.e. 3e14 W/cm^2 (starting values entered here were found empircally by looking at which numbers the code ends up with for a given peak intensity)
        init_grid_size          = 16; % Starting point, will increase with iterations in case it's needed (saddle point not found)
        init_number_points_grid = 100; % Number of points will be increased in case it's needed (if we're in tunneling and that the Coulomb minimum is above the barrier)
    else % "high" intensity
        init_grid_size          = 16; % Starting point, will increase with iterations in case it's needed (saddle point not found)
        init_number_points_grid = 50; % Number of points will be increased in case it's needed (if we're in tunneling and that the Coulomb minimum is above the barrier)
    end
    max_number_points_grid  = 1000; % If the number of points of the grid is above this, something went probably wrong in the code: it exits and returns only NaNs for the given electron
    max_grid_size           = 256;
else
    init_grid_size          = grid_size;
    init_number_points_grid = number_points_grid;
    max_number_points_grid  = grid_size;
    max_grid_size           = number_points_grid;
end

%% "Normal" initialization
conv_elec_field_au      = 5.1422e11;
zci                     = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); % Function to get zero-crossing point
mat_results             = zeros(numel(instants_ionization), 17);
if list_Ip(1) > 0
    list_Ip = - list_Ip; % Potential
end

%% Initialization for coder
BSI                         = 1; % Initialization for coder
x                           = 0; % Initialization for coder
y                           = 0; % Initialization for coder
X                           = zeros(2, 2); % Initialization for coder
Y                           = zeros(2, 2); % Initialization for coder
potential                   = zeros(2, 2); % Initialization for coder
x_saddle                    = NaN; % Initialization for coder
y_saddle                    = NaN; % Initialization for coder
r_saddle                    = NaN; % Initialization for coder
ind_saddle                  = 1; % Initialization for coder
indices_cut_x               = 1; % Initialization for coder
indices_cut_y               = 1; % Initialization for coder
laserField_x                = 1; % Initialization for coder
laserField_y                = 1; % Initialization for coder
multi_x                     = 1; % Initialization for coder
multi_y                     = 1; % Initialization for coder
height_saddle               = NaN; % Initialization for coder
diff_saddle_Ip              = NaN; % Initialization for coder
radiusEntryTunnel           = NaN; % Initialization for coder
radiusExitTunnel            = NaN; % Initialization for coder
xEntryTunnel                = NaN; % Initialization for coder
yEntryTunnel                = NaN; % Initialization for coder
xExitTunnel                 = NaN; % Initialization for coder
yExitTunnel                 = NaN; % Initialization for coder
velocity_after_ionization   = NaN; % Initialization for coder
velocity_after_ionization_x = NaN; % Initialization for coder
velocity_after_ionization_y = NaN; % Initialization for coder
potentialEntryTunnel        = NaN; % Initialization for coder
potentialExitTunnel         = NaN; % Initialization for coder
r                           = NaN; % Initialization for coder
cutProfile                  = 1; % Initialization for coder
indicesTunnelPoints         = 1; % Initialization for coder

%% Going through each ion
for i = 1 : numel(list_Ip)
    codeWentWrong               = false;   
    size_grid                   = init_grid_size; % To be reset for each ion
    number_points_grid          = init_number_points_grid; % To be reset for each ion
    tunnelExitPointNotFound     = true; % Starts by assuming that the tunnel exit point has not been found (for BSI, this is just by-passed, only saddle point is required to exit loop)
    while tunnelExitPointNotFound
        currentIonSaddlePointNotFound  = true;
        while currentIonSaddlePointNotFound % We start from a small grid to maximize the precision in determing the position of the saddle point, and if don't find it, then we increase the grid size
            x                           = linspace(-size_grid, size_grid, number_points_grid); % number of points should be even to prevent zeros from appearing
            y                           = x; % x = y, otherwise, funky problems will appear only sometimes ! (e.g. improfile may for some values return different vector sizes ...)
            [X, Y]                      = meshgrid(x, y);
            indexIonization             = findIndex(t, instants_ionization(i));
            instantaneousAverageCharge  = ionization_level_vs_time(indexIonization);
            laserField_x                = Ex(indexIonization) / conv_elec_field_au;
            laserField_y                = Ey(indexIonization) / conv_elec_field_au;
            laserPotential              = laserField_x .* X + laserField_y .* Y;
            Z_c                         = instantaneousAverageCharge + 1;
            ionPotential                = -Z_c./ (soft_core_param + X.^2 + Y.^2).^(0.5); % https://journals.aps.org/pra/abstract/10.1103/PhysRevA.98.043407
            potential                   = laserPotential + ionPotential;
            [~, ~, nz]                  = surfnorm(X, Y, potential); % normal vectors
            [~, ind_saddle]             = max(reshape(nz, [1 numel(x)*numel(y)])); % The usual x(:) doesn't work with the coder

            if ind_saddle == 1 || ... % Saddle point not found: one needs to increase the size of the map and restart the calculation
                    (abs(X(ind_saddle)) > x(end) - 3) || ... % Avoid edges
                    (abs(Y(ind_saddle)) > y(end) - 3) || ...  % Avoid edges
                    sqrt(X(ind_saddle)^2 + Y(ind_saddle)^2) < 2 * (x(2)-x(1)) % Coulomb center is a saddle point, but not the one we want (2 is just safety, we assume here that x and y are centered around 0)
                currentIonSaddlePointNotFound   = true;
                size_grid                       = 2 * size_grid; % Increases grid size to ensure the saddle point is somewhere within the grid
            else % All good :)
                currentIonSaddlePointNotFound   = false;
                x_saddle                        = X(ind_saddle);
                y_saddle                        = Y(ind_saddle);
                height_saddle                   = potential(ind_saddle);
                diff_saddle_Ip                  = height_saddle - list_Ip(i);

                if diff_saddle_Ip < 0
                    BSI                     = 1;
                    tunnelExitPointNotFound = false; % No tunnel exit point to find: we can exit the main loop!
                else
                    BSI = 0;
                end
            end
            if size_grid > max_grid_size % This should never happen (most cases require a grid size between 4 and 16), but in case something goes wrong, I add this as a safety
                codeWentWrong = true;
                break;
            end
        end % While loop to find saddle point

        if ~codeWentWrong
            if ~BSI % We're in tunnel, let's find the tunnel points
                valuePotentialBelowSaddlePoint = potential(sqrt(X.^2 + Y.^2) < sqrt(x_saddle^2 + y_saddle^2));
                if min(valuePotentialBelowSaddlePoint) > list_Ip(i) % We have so little points that the minimum of the soft-Coulomb is above the barrier. Good luck to find a tunnel point there ... We need to increase the number of points
                    number_points_grid = 2 * number_points_grid; % Increasing RESOLUTION, not grid size
                    if number_points_grid > max_number_points_grid % If have to add so many points, that means that the intensity is really low, and the probability for tunnel ionization becomes really low. This should happen only for intensities ~1e13 W/cm^2 in krypton.
                        codeWentWrong = true; % In principle, we continue to increase the number of points, we will find the tunnel points, but this will be extremely computationnaly expensive, I prefer "losing" the electron
                        break;
                    end
                else % We have enough sampling points to ensure that the Coulomb minimum is below Ip
                    theta = abs(atand(x_saddle / y_saddle));
                    if x_saddle <= 0 && y_saddle <= 0 && theta <= 45 % Bottom left GOOD FOR SURE
                        indices_cut_x = [findIndex(x, -tand(theta) * max(x)), findIndex(x, tand(theta) * max(x))];
                        indices_cut_y = [1 numel(y)];
                    elseif x_saddle <= 0 && y_saddle <= 0 && theta >= 45 % Bottom left GOOD FOR SURE
                        indices_cut_x = [1 numel(x)];
                        indices_cut_y = [findIndex(y, -tand(90-theta) * max(x)) findIndex(y, tand(90-theta) * max(x))];
                    elseif x_saddle <= 0 && y_saddle >= 0 && theta <= 45 % Top left GOOD FOR SURE
                        indices_cut_x = [findIndex(x, tand(theta) * max(x)), findIndex(x, -tand(theta) * max(x))];
                        indices_cut_y = [1, numel(y)];
                    elseif x_saddle <= 0 && y_saddle >= 0 && theta >= 45 % Top left GOOD FOR SURE
                        indices_cut_x = [1, numel(x)];
                        indices_cut_y = [findIndex(y, tand(90-theta) * max(x)) findIndex(y, -tand(90-theta) * max(x))];
                    elseif x_saddle >= 0 && y_saddle >= 0 && theta <= 45 % Top right GOOD FOR SURE
                        indices_cut_x = [findIndex(x, -tand(theta) * max(x)), findIndex(x, tand(theta) * max(x))];
                        indices_cut_y = [1 numel(y)];
                    elseif x_saddle >= 0 && y_saddle >= 0 && theta >= 45 % Top right GOOD FOR SURE
                        indices_cut_x = [1 numel(x)];
                        indices_cut_y = [findIndex(y, -tand(90-theta) * max(x)) findIndex(y, tand(90-theta) * max(x))];
                    elseif x_saddle >= 0 && y_saddle <= 0 && theta <= 45 % Bottom right GOOD FOR SURE
                        indices_cut_x = [findIndex(x, tand(theta) * max(x)), findIndex(x, -tand(theta) * max(x))];
                        indices_cut_y = [1 numel(y)];
                    elseif x_saddle >= 0 && y_saddle <= 0 && theta >= 45 % Bottom right GOOD FOR SURE
                        indices_cut_x = [1 numel(x)];
                        indices_cut_y = [findIndex(y, tand(90-theta) * max(x)) findIndex(y, -tand(90-theta) * max(x))];
                    else
                        error(['I have x_saddle = ', sprintf('%0.5g', x_saddle), ', y_saddle=', sprintf('%0.5g', y_saddle), ' and theta = ', sprintf('%0.5g', theta), '. Something is wrong.'])
                    end

                    %% Get cut profile (improfile, griddedinterpolant and interp2 don't work for mex generation, so I implemented a "home-made" version of it)
                    if diff(indices_cut_x) > diff(indices_cut_y)
                        K               = diff(indices_cut_y) / diff(indices_cut_x);
                        C               = indices_cut_y(1) - K * indices_cut_x(1);
                        xint            = indices_cut_x(1) : indices_cut_x(2);
                        yTarget         = xint * K + C;
                        cutProfile      = zeros(1, numel(xint)); % I interpolate myself since interp2 and griddedinterpolant are not mex compatible ...
                        for k = 1 : numel(xint)
                            y_pix_1                 = floor(yTarget(k));
                            y_pix_2                 = ceil(yTarget(k));
                            dist_pix_1_to_target    = abs(y_pix_1 - yTarget(k));
                            dist_pix_2_to_target    = abs(y_pix_2 - yTarget(k));
                            if dist_pix_1_to_target < 1e-10 || dist_pix_2_to_target < 1e-10
                                cutProfile (k)         = potential(y_pix_2, xint(k)) - list_Ip(i);
                            else
                                cutProfile(k)          = potential(y_pix_1, xint(k)) * (1 - dist_pix_1_to_target) + ...
                                                            potential(y_pix_2, xint(k)) * (1 - dist_pix_2_to_target) - list_Ip(i);
                            end
                        end
                    elseif diff(indices_cut_x) < diff(indices_cut_y)
                        K               = diff(indices_cut_x) / diff(indices_cut_y);
                        C               = indices_cut_x(1) - K * indices_cut_y(1);
                        yint            = indices_cut_y(1) : indices_cut_y(2);
                        xTarget         = yint * K + C;
                        cutProfile      = zeros(1, numel(yint)); % I interpolate myself since interp2 and griddedinterpolant are not mex compatible ...
                        for k = 1 : numel(yint)
                            x_pix_1                 = floor(xTarget(k));
                            x_pix_2                 = ceil(xTarget(k));
                            dist_pix_1_to_target    = abs(x_pix_1 - xTarget(k));
                            dist_pix_2_to_target    = abs(x_pix_2 - xTarget(k));
                            if dist_pix_1_to_target < 1e-10 || dist_pix_2_to_target < 1e-10
                                cutProfile(k)           = potential(yint(k), x_pix_2) - list_Ip(i);
                            else
                                cutProfile(k)           = potential(yint(k), x_pix_1) * (1 - dist_pix_1_to_target) + ...
                                                            potential(yint(k), x_pix_2) * (1 - dist_pix_2_to_target) - list_Ip(i);
                            end
                        end
                    else % Pure diagonal
                        cutProfile = diag(potential - list_Ip(i));
                    end

                    if theta < 45
                        r = x / cosd(theta); % Scaling needed because of diagonal appearing in case of elliptical polarization
                    else
                        r = x / sind(theta); % Scaling needed because of diagonal appearing in case of elliptical polarization
                    end
                    indicesTunnelPoints         = zci(cutProfile);
                    if numel(indicesTunnelPoints) < 4 % With tunnnel, we should have 4 points here: 1 entry, 1 exit, 1 entry on wrong side of barrier, and 1 end of array. If less numbers, no exit for sure !
                        tunnelExitPointNotFound = true;
                        size_grid               = size_grid * 2; % Increases grid size to ensure that we will capture the electron exit position
                    else % We have enough zero crossing-points, let's check if they are meaningful
                        indicesTunnelPoints(indicesTunnelPoints == 1)                   = []; % Removes extreme points which are meaningless
                        indicesTunnelPoints(indicesTunnelPoints == numel(cutProfile))   = []; % Removes extreme points which are meaningless
                        if numel(indicesTunnelPoints) > 3 % 3 points are still possible at this point: entry tunnel, entry tunnel wrong side, and exit tunnel. More means something is wrong
                            error('Too many points in the zero-crossing point detection function');
                        elseif numel(indicesTunnelPoints) < 3 % Electron entered the barrier, but didn't exit ! We need to increase the grid size to get the exit point
                            tunnelExitPointNotFound = true;
                            size_grid               = size_grid * 2; % Increases grid size to ensure that we will capture the electron exit position
                        else % All good: electron entered and exited barrier!
                            if numel(nonzeros(indicesTunnelPoints > numel(cutProfile)/2)) > 1
                                indicesTunnelPoints(1)  = []; % Removes the point from the wrong side of the barrier
                            else
                                indicesTunnelPoints(3)  = []; % Removes the point from the wrong side of the barrier
                            end

                            if x_saddle <= 0 && y_saddle <= 0 && theta <= 45 % Bottom left GOOD FOR SURE
                                multi_x = 1;
                                multi_y = 1;
                            elseif x_saddle <= 0 && y_saddle <= 0 && theta >= 45 % Bottom left GOOD FOR SURE
                                multi_x = 1;
                                multi_y = 1;
                            elseif x_saddle <= 0 && y_saddle >= 0 && theta <= 45 % Top left GOOD FOR SURE
                                multi_x = -1;
                                multi_y = 1;
                            elseif x_saddle <= 0 && y_saddle >= 0 && theta >= 45 % Top left GOOD FOR SURE
                                multi_x = 1;
                                multi_y = -1;
                            elseif x_saddle >= 0 && y_saddle >= 0 && theta <= 45 % Top right GOOD FOR SURE
                                multi_x = 1;
                                multi_y = 1;
                            elseif x_saddle >= 0 && y_saddle >= 0 && theta >= 45 % Top right GOOD FOR SURE
                                multi_x = 1;
                                multi_y = 1;
                            elseif x_saddle >= 0 && y_saddle <= 0 && theta <= 45 % Bottom right GOOD FOR SURE
                                multi_x = -1;
                                multi_y = 1;
                            elseif x_saddle >= 0 && y_saddle <= 0 && theta >= 45 % Bottom right GOOD FOR SURE
                                multi_x = 1;
                                multi_y = -1;
                            end % Else is pointless: there is already a warning above in case something happens

                            radiusEntryTunnel       = r(min(indicesTunnelPoints));
                            radiusExitTunnel        = r(max(indicesTunnelPoints));
                            xEntryTunnel            = multi_x * radiusEntryTunnel * sind(theta);
                            yEntryTunnel            = multi_y * radiusEntryTunnel * cosd(theta);
                            xExitTunnel             = multi_x * radiusExitTunnel * sind(theta);
                            yExitTunnel             = multi_y * radiusExitTunnel * cosd(theta);
                            potentialEntryTunnel    = cutProfile(indicesTunnelPoints(1)) + list_Ip(i);
                            potentialExitTunnel     = cutProfile(indicesTunnelPoints(2)) + list_Ip(i);
                            tunnelExitPointNotFound = false; % Entry and exit points are found!
                            if abs(xEntryTunnel) > abs(xExitTunnel) % Indices in wrong order
                                tempo_storage = xExitTunnel;
                                xExitTunnel     = xEntryTunnel;
                                xEntryTunnel    = tempo_storage;
                            end
                            if abs(yEntryTunnel) > abs(yExitTunnel) % Indices in wrong order
                                tempo_storage = yExitTunnel;
                                yExitTunnel     = yEntryTunnel;
                                yEntryTunnel    = tempo_storage;
                            end
                        end
                    end
                end
            end % If we are in tunnel mode
        end % Did code go wrong in the saddle point detection section?
        if codeWentWrong
            break;
        end
        if ~isnan(BSI)
            if BSI % As explained in supplementary of https://www.nature.com/articles/nphys2125, we have to consider the difference of potential between saddle point and ionization potential as an initial velocity, right at the ionization time
                velocity_after_ionization   = sqrt(2 * abs(diff_saddle_Ip)); % delta_E = 0.5 * m_elec * speed^2, so in atomic units: speed = sqrt(2 * delta_E)
                velocity_after_ionization_x = - velocity_after_ionization * laserField_x / sqrt(laserField_x^2 + laserField_y^2);
                velocity_after_ionization_y = - velocity_after_ionization * laserField_y / sqrt(laserField_x^2 + laserField_y^2);
            else % Tunnel: no velocity at tunnel exit
                velocity_after_ionization_x = 0;
                velocity_after_ionization_y = 0;
            end
        end
    end % While tunnel exit point not found

    if codeWentWrong
        x_saddle                    = NaN;
        y_saddle                    = NaN;
        height_saddle               = NaN;
        diff_saddle_Ip              = NaN;
        BSI                         = NaN;
        velocity_after_ionization   = NaN;
        velocity_after_ionization_x = NaN;
        velocity_after_ionization_y = NaN;
        radiusEntryTunnel           = NaN;
        radiusExitTunnel            = NaN;
        xEntryTunnel                = NaN;
        yEntryTunnel                = NaN;
        xExitTunnel                 = NaN;
        yExitTunnel                 = NaN;
        potentialEntryTunnel        = NaN;
        potentialExitTunnel         = NaN;
    end

    %% From here: only plots
%     if plot_result
%         if BSI
%             multi_factor_display        = 3;
%             potential                   = (potential - height_saddle) * multi_factor_display;
%             potential(potential > 13)   = NaN;
%             potential(potential < -6)   = NaN;
%             list_Ip(i)                  = 5 / multi_factor_display;
%         else
%             multi_factor_display        = 15;
%             potential                   = (potential - potentialEntryTunnel) * multi_factor_display - 1;
%             potential(potential > 15)   = NaN;
%             potential(potential < -6)   = NaN;
%             newHeightEntry              = - 1;
%             newHeightExit               = - 1;
%             list_Ip(i)                  = newHeightEntry / multi_factor_display;
%         end
%         disp(['Showing result from code ', mfilename, ' line 344'])
%         size_grid
%         number_points_grid
%         x_saddle
%         y_saddle
%         height_saddle
%         diff_saddle_Ip
%         BSI
%         if ~isnan(BSI)
%             if ~BSI
%                 radiusEntryTunnel
%                 radiusExitTunnel
%                 xEntryTunnel
%                 yEntryTunnel
%                 xExitTunnel
%                 yExitTunnel
%                 potentialEntryTunnel
%                 potentialExitTunnel
%             end
%         end
%         clf
%         if BSI
%             if ~isnan(velocity_after_ionization) % We're in BSI and we want to show the speed of the ejected electron 
%                 t = quiver3D([x_saddle, y_saddle, 0], [0, 0, list_Ip(i) * multi_factor_display]);
%                 hold on
%                 t = quiver3D([x_saddle, y_saddle, list_Ip(i) * multi_factor_display], [5*velocity_after_ionization_x, 5*velocity_after_ionization_y, 0]);
%             end
%         else
%             t = quiver3D([xEntryTunnel + 1, yEntryTunnel , -1], [xExitTunnel - xEntryTunnel - 2, yExitTunnel  - yEntryTunnel , 0]);
%         end
%         lighting phong
%         camlight head
%         hold on
%         grid on
%         drawnow
% 
%         if ~isnan(x_saddle)
%             if ~BSI % Tunnel: We need a cut profile to see entry and exit points
%                 BSI_or_tunnel = 'tunnel';
%             else
%                 BSI_or_tunnel = 'BSI';
%             end
%             p1 = surf(X, Y, potential);
%             set(p1,'facealpha',0.8)
%             if ~BSI
%                 plot3([xEntryTunnel, xExitTunnel], ...
%                       [yEntryTunnel, yExitTunnel], ...
%                       [newHeightEntry, newHeightExit], ...
%                       '.', 'MarkerSize', 50, 'color', 'red')
%             end
%             p2 = surf(X, Y, list_Ip(i) * ones(size(potential)) * multi_factor_display);
%             set(p2,'facealpha',0.4)
%             
%             hold off
%             if BSI
%                 view(-25, 10)
%             else
%                 view(7, 6)
%             end
%             shading flat
%             colormap(jet)
%             axis square
%             if ~BSI % Tunnel: We need a cut profile to see entry and exit points
% %                     plot(r, cutProfile + list_Ip(i))
% %                     hold on
% %                     plot(r(indicesTunnelPoints), [potentialEntryTunnel, potentialExitTunnel], '*')
% %                     hold off
%             end
% %             title(['Ion number ', sprintf('%1.1f', i), '(', BSI_or_tunnel, ')'])
%         else
%             text(1, 1, 'Code went wrong.')
%         end
%         drawnow
%     end

    %% Let's store the results of the current electron!
    mat_results(i, :)   = [ BSI, ...
                            x_saddle, ...
                            y_saddle, ...
                            sqrt(x_saddle^2 + y_saddle^2), ...
                            height_saddle, ...
                            diff_saddle_Ip, ...
                            velocity_after_ionization, ...
                            velocity_after_ionization_x, ...
                            velocity_after_ionization_y, ...
                            radiusEntryTunnel, ...
                            radiusExitTunnel, ...
                            xEntryTunnel, ...
                            yEntryTunnel, ...
                            xExitTunnel, ...
                            yExitTunnel, ...
                            potentialEntryTunnel, ...
                            potentialExitTunnel];
end % For loop current ion
end % Function