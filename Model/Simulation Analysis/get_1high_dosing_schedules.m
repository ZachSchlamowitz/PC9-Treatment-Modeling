%% get_1high_dosing_schedules.m
% This function outputs all possible dosing schedules of desired length 
% where only one day gets a high dose and all the others get low doses

function treatments = get_1high_dosing_schedules(n, H, L)
        treatments = L*(ones(n)-eye(n)) + H*eye(n);
end
