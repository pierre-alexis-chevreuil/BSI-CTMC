function output = functionGenerator(x, FWHM, center, type, varargin)
%% Function to generate some classical functions (Gaussian, sech^2, ...)
% 
% Inputs : - x : list of numbers for x axis
%          - FWHM : FWHM in x-axis units
%          - center : center of the function in x-axis units
%          - type : string in the following possibilites :
%               - 'Gaussian'
%               - 'Supergaussian'
%               - 'Sech2'
%               - 'Square'
%          - additional parameter : order supergaussian, please look at
%               the example below 
%
% Example :
% 
%     x = 0 : 0.001 : 4;
%     y1 = functionGenerator(x, 1.23, 2, 'Gaussian');
%     y2 = functionGenerator(x, 1.23, 2, 'Supergaussian', 'order', 10);
%     y3 = functionGenerator(x, 1.23, 2, 'Sech2');
%     y4 = functionGenerator(x, 1.23, 2, 'Square');
%     plot(x,y1,x,y2,x,y3,x,y4,'linewidth',3);
%     hold on
%     plot([0 4], [0.5 0.5], 'k--','linewidth',3)
%     legend({'Gaussian', 'Hypergaussian', 'Sech^2', 'Square'})
%     makeItNice
%
% Date : 07.10.2019
%
% Author : Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Decompose varargin
% defaultOrderSuperGaussian   = '10';
% p                           = inputParser; % Creates the parser
% addParameter(p, 'order', defaultOrderSuperGaussian, @(x) isnumeric(x) && isscalar(x) && (x >= 0));
% parse(p, varargin{:});
% superGaussianOrder = p.Results.order;

switch type
    case 'Gaussian'
        output = exp(- (x - center).^2 / ((FWHM)/(2*sqrt(log(2))))^2 );
        
%     case 'Supergaussian' % https://forum.image.sc/t/fwhm-application-to-linear-structures/11258/3
%         sigma = FWHM / (2* sqrt(2) * (log(2)) ^ (1/(2*superGaussianOrder)) );
%         output = exp(- ((x - center).^2 / (2 * sigma^2)) .^ superGaussianOrder);
        
    case 'Sech2' % 10.1364/OL.20.001160
        output = 1 ./ cosh(2*acosh(sqrt(2)) * (x - center) / FWHM) .^ 2;
        
    case 'Square'
        [~, index1] = find(x > center-FWHM/2, 1, 'first');
        [~, index2] = find(x > center+FWHM/2, 1, 'first');
        output = [zeros(1, index1) , ones(1, index2 - index1) , zeros(1, length(x)-index2)];
        
    otherwise
        disp(['You entered an invalid possibility in the ',mfilename]);
        disp('code. This should be one of the following possibilities :');
        disp('"Gaussian", "Supergaussian", "Sech" or "Square".');
end

end