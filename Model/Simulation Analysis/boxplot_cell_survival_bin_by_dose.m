%% boxplot_cell_survival_bin_by_dose.m
% This function plots boxplots of remaining cell counts from
% HML_dose_simulation.m trajectories for difference dosing groups

% Author: Zach Schlamowitz
% Date: 11/29/21 

function boxplot_cell_survival_bin_by_dose(trajectories)
    
    % Create matrices for plotting; one for data, one for binning / groups
    surviving_cell_counts = cell2mat(trajectories(2:end,4));
    total_dose = round(cell2mat(trajectories(2:end,6)),1);
    
    % Identify how many samples fall in each group UNFINISHED
    groups = zeros(size(total_dose,1),1);
    entry = 1;
    for i=1:size(total_dose,1)-1
        if total_dose(i,1) ~= total_dose(i+1,1)
            entry = entry + 1;
        end
        groups(i,1) = entry;
    end
    
    
    group_counts = [];
    count = 0;
    for i=1:size(groups,1)-1
        count = count + 1;
        group_counts(i,1) = count;
        if groups(i,1) ~= groups(i+1,1)
            count = 0;
        end
    end
    
    if groups(end-1,1) ~= groups(end,1)
        group_counts(end+1,1) = 1;
    else
        group_counts(end+1,1) = group_counts(end,1);
    end
    
    data_with_group_counts = [total_dose group_counts];
    disp(data_with_group_counts)
    
    % Plot
    figure('Name','Sim Results Boxplots');
    boxplot(surviving_cell_counts, total_dose);
    title('Cell Counts Remaining After 7 Day Dose Schedules (Simulated)')
    xlabel('Total Dose Administered')
    ylabel('Remaining Cell Count, uM')
%     axis([0.5 2.5 200 400])

end