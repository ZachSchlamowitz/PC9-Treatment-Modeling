%% compute_AUC_onehigh_sims.m

one_high_trajectories = get_1high_dosing_schedules(5, 1.150, 0.2);
AUCs = [];
for i=1:size(one_high_trajectories,1)
    AUCs(i,1) = get_AUC(one_high_trajectories(i,:));
end
