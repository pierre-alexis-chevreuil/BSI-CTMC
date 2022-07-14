function varargout = FWHM(x, y, varargin)
%% Function to compute the FWHM of any list.
% The function will return :
%       - The FWHM in the unit of the x list,
%       - The FWHM in x indices unit
%       - The index of the closest index before
%           the y list reaches 50%
%       - The index of the closest index before
%           the y list goes below 50%
%       - The x value before the y list reaches 50%
%       - The x value before the y list goes below 50%
%
% To get a better precision on the FWHM calculation,
% the code interpolates the data around the points of interest
%
% Inputs : - x : list of numbers
%          - y : list of numbers (same length as x)
%          - additional parameters : 
%               - numberPointsInterpolation : number of points used in the
%                   interpolation to get a more precise value. The bigger,
%                   the better of course.
%               - level : if what is wanted is not the width at half
%                   maximum but 1/e2 for example, you can provide a
%                   different level of evaluation (0.135 for 1/e2 for 
%                   example)
%               - bypassInterp : boolean to bypass the interpolation
% 
% Example of use :
%
%   x       = -10:0.1:10;
%   y       = functionGenerator(x, pi, 1, 'Gaussian');
%   [FWHM_in_x_units, FWHM_in_indices] = FWHM(x,y);
%   FWHM_precise = FWHM(x,y,'numberPointsInterpolation',1e5);
%   disp(['Simple precision : pi = ',num2str(FWHM_in_x_units)]);
%   disp(['Improved precision : pi = ',num2str(FWHM_precise)]);
%
% Date : 17.07.2020 (added possibility to change the level of evaluation)
% 06.01.2021: Added the possibility to bypass the interpolation (useful for
% log plots)
%
% Author : Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Parser
defaultNumberPointsInterpolation    = 100;
defaultLevel                        = 0.5;
defaultBypassInterp                 = false;
p                                   = inputParser; % Creates the parser
addParameter(p, 'numberPointsInterpolation',    defaultNumberPointsInterpolation,   @(x) isnumeric(x) && isscalar(x) && (x >= 0));
addParameter(p, 'level',                        defaultLevel,                       @(x) isnumeric(x) && isscalar(x) && (x >= 0) && (x <= 1));
addParameter(p, 'bypassInterp',                 defaultBypassInterp,                @islogical);
parse(p, varargin{:});
numberPointsInterpolation   = p.Results.numberPointsInterpolation;
level                       = p.Results.level;
bypassInterp                = p.Results.bypassInterp;

%% Classical calculation
y               = y/max(y); % Normalization
index1          = find(y > level, 1, 'first'); % Classical search for index
index2          = find(y > level, 1, 'last'); % Classical search for index
FWHM_in_indices = index2 - index1;

%% Spline interpolation around the 2 points of interest for improved accuracy
if ~bypassInterp
    interpolation_function_1    = @(z) interp1(x(index1-3 : index1+2) , y(index1-3 : index1+2) , z, 'spline');
    interpolation_function_2    = @(z) interp1(x(index2-2 : index2+3) , y(index2-2 : index2+3) , z, 'spline');
    x_interpolated_1            = linspace(x(index1-3) , x(index1+2), numberPointsInterpolation); % Replaces the 5 points around the points of 
    x_interpolated_2            = linspace(x(index2-2) , x(index2+3), numberPointsInterpolation); % interest by 100 points for better accuracy
    y_interpolated_1            = interpolation_function_1(x_interpolated_1);
    y_interpolated_2            = interpolation_function_2(x_interpolated_2);
    [~, index1_interpolated]    = find(y_interpolated_1 > level, 1, 'first');
    [~, index2_interpolated]    = find(y_interpolated_2 > level, 1, 'last');
    index1_in_x_units           = x(index1-3) + (x_interpolated_1(2) - x_interpolated_1(1) ) * index1_interpolated;
    index2_in_x_units           = x(index2-2) + (x_interpolated_2(2) - x_interpolated_2(1) ) * index2_interpolated;
    FWHM_in_x_units             = index2_in_x_units - index1_in_x_units;
else
    index1_in_x_units   = x(index1);
    index2_in_x_units   = x(index2);
    FWHM_in_x_units     = index2_in_x_units - index1_in_x_units;
end

%% Output
varargout{1} = FWHM_in_x_units;
varargout{2} = FWHM_in_indices;
varargout{3} = index1;
varargout{4} = index2;
varargout{5} = index1_in_x_units;
varargout{6} = index2_in_x_units;
end