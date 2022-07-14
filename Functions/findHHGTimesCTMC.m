function output_mat = findHHGTimesCTMC(t, x, vx, y, vy, z, vz, radius_recombination_sphere, Ip_au_current_ion, lambda, BSI) %#codegen
%% Code to get various parameters about a specific trajectory. It will tell
% if it recombines, and if yes, at what time, with which energy ...
%
% Outputs:  - 1: max distance after 1 cycle (meters)
%           - 2: max distance after 2 cycles (meters)
%           - 3: max distance after 3 cycles (meters)
%           - 4: max energy after 1 cycle (eV), accounts for Ip of ion
%           - 5: max energy after 2 cycles (eV), accounts for Ip of ion
%           - 6: max energy after 3 cycles (eV), accounts for Ip of ion
%           - 7: energy after 1 cycle (eV), accounts for Ip of ion
%           - 8: energy after 2 cycles (eV), accounts for Ip of ion
%           - 9: energy after 3 cycles (eV), accounts for Ip of ion
%           - 10: energy at end of the pulse (eV), accounts for Ip of ion
%           - 11: time of recombination (seconds), taking 0 time at emission
%           - 12: energy of recombination (eV), accounts for Ip of ion
%           - 13: Ip of ion (in a.u.)
%           - 14: BSI or not BSI (bool)
%
% Date: 25.04.2022 (16.05.2022 : modif to account for z and v_z)
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

t                                               = t - t(1); % Removes offset of starting time (start time doesn't matter)
time_after_emission_recombination_not_possible  = 0.33 / (2 * pi) * lambda / 299792458; % If electron recombines before this time, it's before the moment it can recombine in principle, so it's artefact. 0.33 rad is phase of 
distanceToZero                                  = sqrt(x.^2 + y.^2 + z.^2);
speed                                           = sqrt(vx.^2 + vy.^2 + vz.^2);
indicesOutsideSphere                            = distanceToZero > radius_recombination_sphere; % Equivalent of find, but expected to be faster
indicesHHGPossible                              = t > time_after_emission_recombination_not_possible;
one_cycle                                       = lambda / 299792458;
allIndices                                      = (0 : numel(t) - 1)';
indicesOneCycle                                 = nonzeros(allIndices .* (t <= one_cycle));
indicesTwoCycles                                = nonzeros(allIndices .*  ((t > one_cycle) .* (t < 2*one_cycle)));
indicesThreeCycles                              = nonzeros(allIndices .* ((t > 2*one_cycle) .* (t < 3*one_cycle)));
indicesHHG                                      = nonzeros(allIndices .* ~indicesOutsideSphere .* indicesHHGPossible);  
m_e                                             = 9.10938356e-31;
conversion_J_to_eV                              = 1.602176634e-19;
Ip_current_ion_eV                               = Ip_au_current_ion * 27.2113962; % From atomic units to eV

if any(isnan(x)) || any(distanceToZero == 0)
    output_mat = [  NaN, ... % 1: max distance after 1 cycle (meters)
                    NaN, ... % 2: max distance after 2 cycles (meters)
                    NaN, ... % 3: max distance after 3 cycles (meters)
                    NaN, ... % 4: max energy after 1 cycle (eV), accounts for Ip of ion
                    NaN, ... % 5: max energy after 2 cycles (eV), accounts for Ip of ion
                    NaN, ... % 6: max energy after 3 cycles (eV), accounts for Ip of ion
                    NaN, ... % 7: energy after 1 cycle (eV), accounts for Ip of ion
                    NaN, ... % 8: energy after 2 cycles (eV), accounts for Ip of ion
                    NaN, ... % 9: energy after 3 cycles (eV), accounts for Ip of ion
                    NaN, ... % 10: energy at end of the pulse (eV), accounts for Ip of ion
                    NaN, ... % 11: time of recombination (seconds), taking 0 time at emission
                    NaN, ... % 12: energy of recombination (eV), accounts for Ip of ion
                    NaN, ... % 13: origin of ion (number, so 0, 1, 2, ...)
                    NaN]; % 14: BSI or not BSI (bool)
else
    max_distance_one_cycle      = max(distanceToZero(indicesOneCycle));
    max_distance_two_cycles     = max(distanceToZero(indicesTwoCycles));
    max_distance_three_cycles   = max(distanceToZero(indicesThreeCycles));
    max_energy_one_cycle        = 0.5 * m_e * max(speed(indicesOneCycle)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    max_energy_two_cycles       = 0.5 * m_e * max(speed(indicesTwoCycles)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    max_energy_three_cycles     = 0.5 * m_e * max(speed(indicesThreeCycles)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    energy_after_one_cycle      = 0.5 * m_e * speed(indicesOneCycle(end)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    energy_after_two_cycles     = 0.5 * m_e * speed(indicesTwoCycles(end)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    energy_after_three_cycles   = 0.5 * m_e * speed(indicesThreeCycles(end)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    energy_end_pulse            = 0.5 * m_e * speed(end).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    if isempty(indicesHHG)
        time_HHG                = NaN;
        energy_recomb_HHG       = NaN;
    else
        time_HHG                = t(indicesHHG(1));
        energy_recomb_HHG       = 0.5 * m_e * speed(indicesHHG(1)).^2 / conversion_J_to_eV + Ip_current_ion_eV;
    end
    output_mat = [  max_distance_one_cycle, ... % 1: max distance after 1 cycle (meters)
                    max_distance_two_cycles, ... % 2: max distance after 2 cycles (meters)
                    max_distance_three_cycles, ... % 3: max distance after 3 cycles (meters)
                    max_energy_one_cycle, ... % 4: max energy after 1 cycle (eV), accounts for Ip of ion
                    max_energy_two_cycles, ... % 5: max energy after 2 cycles (eV), accounts for Ip of ion
                    max_energy_three_cycles, ... % 6: max energy after 3 cycles (eV), accounts for Ip of ion
                    energy_after_one_cycle, ... % 7: energy after 1 cycle (eV), accounts for Ip of ion
                    energy_after_two_cycles, ... % 8: energy after 2 cycles (eV), accounts for Ip of ion
                    energy_after_three_cycles, ... % 9: energy after 3 cycles (eV), accounts for Ip of ion
                    energy_end_pulse, ... % 10: energy at end of the pulse (eV), accounts for Ip of ion
                    time_HHG, ... % 11: time of recombination (seconds), taking 0 time at emission
                    energy_recomb_HHG, ... % 12: energy of recombination (eV), accounts for Ip of ion
                    Ip_au_current_ion, ... % 13: Ip of ion (in a.u.)
                    BSI]; % 14: BSI or not BSI (bool)
end
end