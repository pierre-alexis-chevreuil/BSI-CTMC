function res_CTMC = get_CTMC_one_HWP_and_QWP_angle_randomized_CEP(params)
%% Same code as get_CTMC_one_HWP_and_QWP_angle, except that this time it
% randomly distributes the CEP (first CTMC calculations have showed a
% "spicky" spectrum, and I think it's due to the fixed CEP.
% 
% Date: 30.04.2022
%
% Author: Pierre-Alexis Chevreuil (chpierre@phys.ethz.ch)

list_CEPs               = params.CEP;
params.numTrajectories  = round(params.numTrajectories / numel(params.CEP));

params.CEP              = list_CEPs(1);
res_CTMC                = get_CTMC_one_HWP_and_QWP_angle(params);

for i = 2 : numel(list_CEPs)
    params.CEP                              = list_CEPs(i);
    res_CTMC.CTMC.CTMC_trajectory_analysis  = [res_CTMC.CTMC.CTMC_trajectory_analysis; get_CTMC_one_HWP_and_QWP_angle(params).CTMC.CTMC_trajectory_analysis];
end
end
