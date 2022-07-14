function [list_Ip, Z_c, l, m] = getIonizationPotential(nameAtom, ionization_level, type_output_Ip)
%% Function to get the ionization potential of noble gases in different
% units and their related properties (charge of the ion, orbital quantum
% number, and magnetic quantum number)
%
% Inputs:   - nameAtom: 'Ar', 'He', ...
%           - ionization_level: 0, 1, 2, ... ground-state is 0 (He = He, 0)
%           - type_output_Ip: 'J', 'eV', 'au'
%
% Outputs:  - Ip_J: Ionization potential in joules, eV, or atomic units,
%               depending on type_output_Ip
%           - Z_c: Charge of the ion (useful for PPT or ADK)
%           - l: Orbital quantum number
%           - m: Magnetic quantum number
%
% Example:
%   
% 
% For calculation of l and m, see on OneNote keV beamline/HHG
% sim/Ionization/Reminders ionization electron levels. Don't forget that
% Hund's rule doesn't apply for ionization calculations (see original ADK
% paper, and review paper from Lai 2017)
%
% Date: 19.04.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

list_Ip = zeros(1, numel(ionization_level));
for i = 1 : numel(ionization_level)
    if strcmp(nameAtom, 'H')
        nameAtomCorrected = 'H';
    else
        if ionization_level(i) == 0
            nameAtomCorrected = nameAtom;
        else
            nameAtomCorrected = [nameAtom, num2str(ionization_level(i))];
        end

        switch nameAtomCorrected
            case 'Ar'
                Ip              = 0.579154843219695; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 1; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Ar1'
                Ip              = 1.01537127301097; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 2; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Ar2'
                Ip              = 1.49698309122411; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 3; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ar3'
                Ip              = 2.18952381429072; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 4; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ar4'
                Ip              = 2.75031826555081; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 5; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ar5'
                Ip              = 3.35484439420275; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 6; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ar6'
                Ip              = 4.5719814994278; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 7; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Ar7'
                Ip              = 5.27193455806579; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 8; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'H'
                Ip              = 0.499733200055859; % I_p (atomic units), Lide2001_Handbook_of_chemistry_and_physics_3rd_edition_9780849304828, page 10-175 (904 in pdf)
                Z_c             = 1; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'He'
                Ip              = 0.903570651554921; % I_p (atomic units), Lide2001_Handbook_of_chemistry_and_physics_3rd_edition_9780849304828, page 10-175 (904 in pdf)
                Z_c             = 1; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'He1'
                Ip              = 1.999816529303915; % I_p (atomic units), Lide2001_Handbook_of_chemistry_and_physics_3rd_edition_9780849304828, page 10-175 (904 in pdf)
                Z_c             = 2; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Kr'
                Ip              = 0.514475824654672; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 1; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Kr1'
                Ip              = 0.895207280837725; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 2; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Kr2'
                Ip              = 1.31702172636037; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 3; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Kr3'
                Ip              = 1.86870234905477; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 4; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Kr4'
                Ip              = 2.37731278191451; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 5; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Kr5'
                Ip              = 2.88445324242495; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 6; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Kr6'
                Ip              = 4.01045206199306; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 7; % Z_c
                l               = 2; % l
                m               = 0; % m
            case 'Kr7'
                Ip              = 4.62313653718364; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 8; % Z_c
                l               = 2; % l
                m               = 0; % m
            case 'Kr8'
                Ip              = 8.56258893470523; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 9; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Kr9'
                Ip              = 9.84881474034765; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 10; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Kr10'
                Ip              = 11.3187870896533; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 11; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Kr11'
                Ip              = 12.8622580564242; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 12; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Kr12'
                Ip              = 14.3689797144624; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 13; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Kr13'
                Ip              = 16.3901916947577; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 14; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Kr14'
                Ip              = 18.0806598964591; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 15; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Kr15'
                Ip              = 19.8446267156259; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 16; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Kr16'
                Ip              = 21.7188414609905; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 17; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Kr17'
                Ip              = 23.5195575888899; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 18; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Mg'
                Ip              = 0.280998125757966; % I_p (atomic units), Lide2001_Handbook_of_chemistry_and_physics_3rd_edition_9780849304828, page 10-175 (904 in pdf)
                Z_c             = 1; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Ne'
                Ip              = 0.792481974886684; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 1; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Ne1'
                Ip              = 1.50536083113589; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 2; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Ne2'
                Ip              = 2.33076243254288; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 3; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ne3'
                Ip              = 3.57166531572533; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 4; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ne4'
                Ip              = 4.63948997956966; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 5; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ne5'
                Ip              = 5.80396532538084; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 6; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Ne6'
                Ip              = 7.61706597032313; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 7; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Ne7'
                Ip              = 8.78664947004814; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 8; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Xe'
                Ip              = 0.445763371009974; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 1; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Xe1'
                Ip              = 0.770816750667134; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 2; % Z_c
                l               = 1; % l
                m               = 0; % m
            case 'Xe2'
                Ip              = 1.14106603614849; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 3; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Xe3'
                Ip              = 1.55082082851743; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 4; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Xe4'
                Ip              = 1.98813760243585; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 5; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Xe5'
                Ip              = 2.45128914039332; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 6; % Z_c
                l               = 1; % l
                m               = 1; % m
            case 'Xe6'
                Ip              = 3.36623667990987; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 7; % Z_c
                l               = 2; % l
                m               = 0; % m
            case 'Xe7'
                Ip              = 3.89461089100603; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 8; % Z_c
                l               = 2; % l
                m               = 0; % m
            case 'Xe8'
                Ip              = 6.60899568247806; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 9; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Xe9'
                Ip              = 7.42336036399338; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 10; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Xe10'
                Ip              = 8.41632668594932; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 11; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Xe11'
                Ip              = 9.37107372682332; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 12; % Z_c
                l               = 2; % l
                m               = 1; % m
            case 'Xe12'
                Ip              = 10.326555753872; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 13; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Xe13'
                Ip              = 11.5392829420491; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 14; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Xe14'
                Ip              = 12.6050128952957; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 15; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Xe15'
                Ip              = 13.7442414660075; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 16; % Z_c
                l               = 2; % l
                m               = 2; % m
            case 'Xe16'
                Ip              = 14.8467207279868; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 17; % Z_c
                l               = 0; % l
                m               = 0; % m
            case 'Xe17'
                Ip              = 15.949199989966; % I_p (atomic units), https://physics.nist.gov/PhysRefData/ASD/ionEnergy.html
                Z_c             = 18; % Z_c
                l               = 0; % l
                m               = 0; % m
            otherwise
                error('%s is not a species known by the ADK_fast code.', nameAtomCorrected);
        end
        if strcmp(type_output_Ip, 'J')
            list_Ip(i) = Ip * 4.3597447222071e-18;
        elseif strcmp(type_output_Ip, 'eV')
            list_Ip(i) = Ip * 27.2113962;
        elseif strcmp(type_output_Ip, 'au')
            list_Ip(i) = Ip;
        else
            error('type_output_Ip should be J, eV or au');
        end
    end
end

end