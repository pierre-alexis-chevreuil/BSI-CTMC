function mat_BSI_vs_ions = getBSIIntensities(t, I_enveloppe_vs_time, gas, max_ion)
%% Code to get the Barrier Suppression Ionizations intensities and the
% times at which such intensities are reached. This is a way to check that
% no significant amount of population is left when the BSI intensnties of
% each individual atom/ion is left, since the Krainov BSI ionization rates
% are inaccurate when the intensity is ~ 10 times the BSI value (this is
% the case for all S-matrix or Keldysh-Faisal-Reiss theories).
%
% Output: mat_BSI_vs_ions:  column 1: ion level
%                           column 2: BSI intensity for the ion
%                           column 3: time at which the BSI intensity is
%                               reached
%
% Date: 17.04.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

mat_BSI_vs_ions                 = zeros(max_ion + 1, 3);
[~, ~, mat_BSI_vs_ions(1, 2)]   = ADK_fast(1e10, gas, 1, 'lin');
mat_BSI_vs_ions(1, 3)           = t(findIndex(I_enveloppe_vs_time, mat_BSI_vs_ions(1, 2)));
for i = 1 : max_ion
    mat_BSI_vs_ions(i+1, 1)         = i;
    [~, ~, mat_BSI_vs_ions(i+1, 2)] = ADK_fast(1e10, [gas, num2str(i)], 1, 'lin'); % linear or circular: BSI intensity is the same
    mat_BSI_vs_ions(i+1, 3)         = t(findIndex(I_enveloppe_vs_time, mat_BSI_vs_ions(i+1, 2)));
end
end