function [] = checkChirpingFFT(time, E)
%% Checks that the pulse duration is not too long compared to the length of
% the length of the time element. This code is used in the waveplates
% simulation.
%
% Date: 02.03.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

if numel(nonzeros(real(E(1, :)))) > 10 % Checks that it's not zero, otherwise code not happy
    y       = real(E(1, :)) ./ max(real(E(1, :)));
elseif numel(nonzeros(real(E(2, :)))) > 10 % Power is in the other polarization
    y       = real(E(1, :)) ./ max(real(E(1, :)));
else
    error('Both polarizations of the E vector have mostly zero values, there must be an error somewhere.')
end
index1      = find(y > 0.5, 1, 'first'); % Classical search for index
index2      = find(y > 0.5, 1, 'last'); % Classical search for index
tau_FWHM    = round((time(index2) - time(index1))*1e15);
tau_FWHM    = tau_FWHM(1); % Because of compiler ... I don't really get why it's needed ...
delta_t     = round((time(end) - time(1))*1e15);
delta_t     = delta_t(1); % Because of compiler ... I don't really get why it's needed ...
if delta_t < 3 * tau_FWHM
    error(['The pulse duration at that part of the system is quite long (', sprintf('%4.0f', tau_FWHM), ...
            ' fs, while the time vector is ', sprintf('%4.2f',delta_t), ' fs. This could lead ' ...
            ' to funky business on the edges because of the FFTs. To prevent this, just make the time axis longer ', ...
            '(or decrease the dispersion).'])
end
end