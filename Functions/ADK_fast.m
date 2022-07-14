function [w_ADK, E_cr_AU, I_BSI_SI] = ADK_fast(E_field_list, nameAtom, BSI_correction, linear_or_circular) %#codegen
% This code is optimized for speed, hence why it doesn't have many
% functions, and that all variables are defined only locally and all the
% things which could be pre-compute have been pre-computed.
%
% Reminder: from I_peak (W/m^2) to E_peak (V/m):
%           E_peak = sqrt(753.460626923542 * I_peak)
%
% Source: https://journals.aps.org/pra/abstract/10.1103/PhysRevA.96.063417
%
% Inputs: - E_field_list: Electric field in SI (V/m), can be one element
%               or a list
%         - nameAtom: 'He', 'Ne', 'Ar', 'Kr', 'Xe', and some of their ions
%               (see code to see which ones)
%         - BSI_correction: ADK following Krainov's extension of ADK
%               (Krainov1997 10.1364/JOSAB.14.000425). Note that I have a
%               doubt on the exponent on the (2F)^x in the Airy function.
%               This is reported as 2/3 in the original publication and
%               another they published the next year (Tunneling and
%               barrier-suppression ionization of atoms and ions in a laser
%               radiation field), but this in opposition with 
%
%
% Output: - w_ADK: Ionizations / second (could be > 1). To get ionization
%               probability : P = 1 - exp(-trapz(list_time, w_PPT_at_fixed_I)),
%               with a new peak intensity at every peak intensity the
%               ionization probability needs to be computed.
%         - E_cr_AU: BSI field strength in atomic units
%         - I_BSI_SI: BSI intensity in W/m^2
%
% Example:
%     tau_FWHM        = 10e-15;
%     time_list       = linspace(-2*tau_FWHM, 2*tau_FWHM, 1000);
%     lambda          = 800e-9;
%     atom            = 'He';
%     I_peak          = 2e15*1e4;
%     f_t             = functionGenerator(time_list, tau_FWHM, 0, 'Gaussian'); % To be scaled later in W/m^2
%     I_field         = I_peak * f_t;
%     CEP             = 0;
%     E_field         = sqrt(753.460626923542 * I_peak * f_t).* cos(2*pi*299792458/lambda*time_list + CEP);
%     E_field_CA      = sqrt(753.460626923542 * I_peak * f_t);
%     E_max           = sqrt(753.460626923542 * I_peak); % Don't take max(E_field)! If CEP ~= 0, this will lead to wrong results!
%     w_YI                = PPT_Yudin_Ivanov(E_max, lambda, atom, time_list, CEP, sqrt(f_t), 'approx');
%     w_PPT               = PPT_fast(E_field, lambda, atom, 0, 'approx');
%     w_PPT_CA            = PPT_fast(E_field_CA, lambda, atom, 1, 'approx');
%     w_ADK               = ADK_fast(E_field, atom, 0);
%     w_ADK_BSI           = ADK_fast(E_field, atom, 1);
%     Ion_prob_YI         = 1 - exp(-cumtrapz(time_list, w_YI));
%     Ion_prob_PPT        = 1 - exp(-cumtrapz(time_list, w_PPT));
%     Ion_prob_PPT_CA     = 1 - exp(-cumtrapz(time_list, w_PPT_CA));
%     Ion_prob_ADK        = 1 - exp(-cumtrapz(time_list, w_ADK));
%     Ion_prob_ADK_BSI    = 1 - exp(-cumtrapz(time_list, w_ADK_BSI));
%     yyaxis left
%       plot(time_list * 1e15, Ion_prob_PPT*100, '-', 'color', [0 0.4470 0.7410])
%       hold on
%       plot(time_list * 1e15, Ion_prob_PPT_CA*100, '-', 'color', [0.9290 0.6940 0.1250])
%       plot(time_list * 1e15, Ion_prob_YI*100, '--', 'color', [0.4940 0.1840 0.5560])
%       plot(time_list * 1e15, Ion_prob_ADK*100, '-', 'color', [0.4660 0.6740 0.1880])
%       plot(time_list * 1e15, Ion_prob_ADK_BSI*100, '-', 'color', [0.3010 0.7450 0.9330])
%       hold off
%       ylabel('Ionization probability (%)')
%       grid on
%     yyaxis right
%       plot(time_list * 1e15, E_field * 1e-10)
%       xlim([time_list(1) time_list(end)]*1e15)
%       ylabel('Electric field (10^{10} V/m)')
%       legend({'PPT', 'PPT CA', 'Yudin-Ivanov', 'ADK', 'ADK BSI corrected', 'Electric field'}, ...
%           'interpreter', 'none', 'location', 'northwest')
%     title(['Ionization probability calculation in ', atom, ' with \lambda=', ...
%       num2str(lambda*1e9), ' nm, I_{peak}=', num2str(I_peak*1e-12*1e-4), ' TW.cm^{-2} and {\tau}_{FWHM}=', ...
%       num2str(tau_FWHM*1e15), ' fs'])
%     makeItNice('screenNumber', 2)
%
% List of modifications:
%   - 06.06.2021: Corrected BSI correction from Krainov 1997 + changed
%       integration by trapz to go 5 x faster + mex generation)
%   - 17.04.2022: Added possibility to compute ionization for circular
%       polarization in BSI and non-BSI (only for s-states in non_BSI !!!)
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

%% Atom properties
switch nameAtom
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
        error('%s is not a species known by the ADK_fast code.', nameAtom);
end

nStar       = Z_c / sqrt(2 * Ip);
F_0         = abs(E_field_list) / 5.14220674763e11; % From V/m to atomic units
E_cr_AU     = (Ip / 0.499733200055859)^2 / (16*Z_c); % Above-barrier field strength in atomic units
I_BSI_SI     = (E_cr_AU* 5.14220826e11).^2  / 753.460626923542;

if nargin < 4
    linear_or_circular = 'lin';
end

if BSI_correction % Krainov1997 10.1364/JOSAB.14.000425, time averaging coming from Lai2017
    w_ADK = zeros(1, length(F_0));
    D     = (4 * exp(1) * Z_c.^3 ./ (F_0 .* nStar.^4)).^nStar; % Equation 8 in Delone1998

    if strcmp(linear_or_circular, 'lin')
        for k = 1 : length(F_0)
            airyFun     = @(x)   x.^2 .* airy(x.^2 + 2 * Ip ./ (2 * F_0(k)).^(2/3) ).^2;
            xx          = linspace(0, 2, 20);
            integ       = trapz(xx, airyFun(xx)); % Can be replaced by integral(airyFun, 0, inf);  but it is much slower (5 times than this with 100 points)
            w_ADK(1, k) = real(4 * sqrt(3) * F_0(k) .* D(k).^2 ./ (pi * nStar .* (2 * F_0(k)).^(1/3)) .* ...
                            integ); % Delone1998 equation 48, equivalent to equation 23 in Krainov1997
        end
    elseif strcmp(linear_or_circular, 'circ') % Delone 10.1070/pu1998v041n05abeh000393
        D           = (4 * exp(1) * Z_c.^3 ./ (F_0 .* nStar.^4)).^nStar; % Equation 8 in Delone1998
        for k = 1 : length(F_0)
            xx          = linspace(-2, 2, 20);
            [X, Y]      = meshgrid(xx, xx);
            I           = airy(X.^2 + Y.^2 + Z_c^2 ./ (nStar.^2 .* (2 * F_0(k)).^(2/3) )).^2;
            integ       = trapz(xx, trapz(xx, I, 2));
            w_ADK(1, k) = real(Z_c .* (2 * F_0(k)).^(1/3) .* D(k).^2 / (2 * pi * nStar.^2) .* ...
                            integ); % equation 14-15 in Krainov1997
        end
    else
        error('The only options possible are ''lin'' or ''circ'' for argument linear_or_circular');
    end

else % Normal ADK
    if strcmp(linear_or_circular, 'lin')
%         D           = (4 * exp(1) * Z_c.^3 ./ (F_0 .* nStar.^4)).^nStar; % Equation 8 in Delone1998
%         w_ADK       = F_0 .* D.^2 ./ (8 * pi * Z_c) .* ...
%                         sqrt(3 * nStar.^3 .* F_0 ./ (pi * Z_c^3)) .* ...
%                         exp(-2 * Z_c^3 ./ (3 * nStar^3 .* F_0)); % Delone1998, equation 7
        lStar       = nStar - 1;
        c_nl        = 2^(2*nStar) / (nStar * gamma(nStar + lStar + 1) * gamma(nStar - lStar));
        f_lm        = (2*l+1) * factorial(l+abs(m)) / (2^abs(m) * factorial(abs(m)) * factorial(l-abs(m)));
        CC          = (2 ./ (F_0 * nStar^3)).^(2 * nStar - abs(m) - 1); % Long range Coulomb correction
        w_ADK       = real(c_nl .* f_lm * Ip .* CC .* exp(- 2 * (2*Ip)^(3/2) ./ (3 * F_0))); % Equation A6 in Lai2017
    elseif strcmp(linear_or_circular, 'circ')
        D           = (4 * exp(1) * Z_c.^3 ./ (F_0 .* nStar.^4)).^nStar; % Equation 8 in Delone1998
        w_ADK       = real(F_0 .* D.^2 ./ (8 * pi * Z_c) .* ...
                            exp(-(2 * Z_c^3 ./ (3 * nStar^3 .* F_0)))); % Delone1998, equation 9
    else
        error('The only options possible are ''lin'' or ''circ'' for argument linear_or_circular');
    end
end
w_ADK = w_ADK / 2.418884326505e-17; % Conversion atomic units to ionizations / second
end