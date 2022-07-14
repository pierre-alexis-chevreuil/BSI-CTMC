function varargout = getDispersionValues(GD_or_phase , lambda_grid , input_GD_or_phase , varargin)
%% getDispersionValues allows you to fit a spectral phase or a spectral
% group delay and get the subsequent dispersion values
%
% Inputs : - GD_or_phase : Array of group delays or spectral phase,
%                       depending on the case (you precise it after).
%          - lambda_grid : Array of wavelength for x axis (in meters).
%          - input_GD_or_phase : Precise if your input was a GD or a
%                       spectral phase. 2 possibilites : 'GD' or
%                       'phase'.
%          - optional parameters :
%            - lambda_central : Wavelength around which the fit will
%                       be done. If nothing is given, the code will
%                       take the mean value of the wavelength array.
%            - order_fit : Order of the fit. The length of the
%                       dispersion array will be equal to it.
%            - displayValues : Boolean to display in the command window the
%                       dispersion values computed (in fs^n)
%            - displayName : Name to display if displayValues is true            
%
% Outputs : - dispersion_values : An array with all the dispersion
%                   values in SI (GD in s, GDD in s^2, TOD in s^3 ...)
%           - lambda_central : only useful if you didn't provide it
%                   in the inputs. It is just the mean of the array
%                   of wavelengths
% 
% Example :
%
%     lambda              = (655 : 1 : 765) * 1e-9;
%     lambda_central      = 730e-9;
%     GD_offset           = 1000      * 1e-15 ^ 1;
%     GDD                 = -2323     * 1e-15 ^ 2;
%     TOD                 = -10365    * 1e-15 ^ 3;
%     FOD                 = -5000     * 1e-15 ^ 4;
%     GD                  = getGDCurve(lambda , lambda_central , [GD_offset GDD TOD FOD]);
%     phase               = GDToPhase(lambda,   GD , 'wavelength');
%     dispersion_values_1 = getDispersionValues(GD ,    lambda , 'GD' ,    lambda_central, 4, ...
%                                               'displayValues', true, 'displayName', 'input is GD');
%     dispersion_values_2 = getDispersionValues(phase , lambda , 'phase' , lambda_central, 4, ...
%                                               'displayValues', true, 'displayName', 'input is phase');
%
% Date : 26.07.2020
%
% Author : Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch) adapted from
% Justinas Pupeikis code

%% Parser
defaultLambdaCentral    = mean(lambda_grid);
defaulOrderFit          = 4;
defaultDisplayValue     = false;
p                       = inputParser; % Creates the parser
addRequired( p, 'GD_or_phase');
addRequired( p, 'lambda_grid');
addRequired( p, 'input_GD_or_phase');

addOptional( p, 'lambda_central',   defaultLambdaCentral,   @isnumeric);
addOptional( p, 'order_fit',        defaulOrderFit,         @isnumeric);

addParameter(p, 'displayValues',    defaultDisplayValue,    @islogical);
addParameter(p, 'displayName',      '',                     @ischar);

parse(p , GD_or_phase , lambda_grid , input_GD_or_phase , varargin{:});

lambda_central      = p.Results.lambda_central;
lambda_grid         = p.Results.lambda_grid;
order_fit           = p.Results.order_fit;
displayValues       = p.Results.displayValues;
displayName         = p.Results.displayName;

if strcmp(input_GD_or_phase , 'GD')
    GD = p.Results.GD_or_phase;
else
    phase = p.Results.GD_or_phase;
end

%% Code
c           = 299792458;
omega       = 2 *pi *c ./ lambda_grid;
omega_0     = 2 *pi *c / lambda_central;
omega_fit   = omega - omega_0;
scale       = 1e15;

if strcmp(input_GD_or_phase , 'GD')
    [p_fit,~]   = polyfit( omega_fit /scale , GD * scale, order_fit); 
    poly        = fliplr(p_fit);
    for i = 1 : length(poly)
        poly(i) = poly(i) * factorial(i - 1);
    end
    for i = 1 : order_fit - 1
        dispersion_values(i) = poly(i+1) / scale^(i + 1);
    end
    [~ , index_omega_0] = min(abs(lambda_grid - lambda_central));
    dispersion_values   = [GD(index_omega_0) dispersion_values];
else
    [p_fit,~]   = polyfit( omega_fit /scale , phase * scale, order_fit + 1);
    poly        = fliplr(p_fit);
    for i = 1 : length(poly)
        poly(i) = poly(i) * factorial(i - 1);
    end
    for i = 1 : order_fit - 0
        dispersion_values(i) = poly(i+1) / scale^(i + 1);
    end
end

if displayValues
    
    disp(pad([' DISPERSION VALUES ', displayName, ' @ ', num2str(lambda_central*1e9, '%4.0f'), ' nm '], 60, 'both', '='))
    for i = 1 : length(dispersion_values)
        disp([num2str(dispersion_values(i) * 1e15^(i) , '%4.0f'),' fs^',num2str(i)])
    end
    disp(pad([' END DISPERSION VALUES ', displayName, ' '], 60, 'both', '='))
end

%% Output parser
varargout{1} = dispersion_values;
varargout{2} = lambda_central;

end