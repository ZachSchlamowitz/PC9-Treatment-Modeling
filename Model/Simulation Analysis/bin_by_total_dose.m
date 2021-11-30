%% bin_by_total_dose.m
% This function takes in the simulated trajectories resulting from
% HML_dose_simulation and bins them by total dose administered.

% Author: Zach Schlamowitz
% Date: 11/29/21

function bins = bin_by_total_dose(trajectories)
    
    % create helper cell array to hold each distinct dose level
    distinct_doses = {};
    k=1;
    for i=1:size(trajectories,1)
        
        cur_dose = trajectories(i,6);
        
        % check if already have that dose noted
        already_have_dose = false;
        for ii=1:size(distinct_doses,1)
            if cellfun(@isequal, cur_dose, distinct_doses(ii,1))
                already_have_dose = true;
                break
            end
        end
        
        if ~already_have_dose
            distinct_doses(k,1) = cur_dose;
            k = k + 1;
        end
        
    end
    
    % sort the distinct doses to increasing order
    
    
    % Initialize bins struct
    bins = struct();
    
    % now that know which distinct doses exist, put all the corresponding
    % treatment schedules and their outcomes together in their bins
    for i=1:size(trajectories,1)
        
        dose = trajectories(i,6);
        
        % get index in distinct doses our current dose has
        j = 0;  
        for ii=1:size(distinct_doses,1)
            if cellfun(@isequal, distinct_doses(ii,1), dose)
                j = ii;
                break
            end
        end
        
        if j >= size(bins,2)
            bins(j).dose = dose;
        end
        
        
        % Set dose of current bin
        if ~isfield(bins(j), 'dose')
            bins(j).dose = dose;
        end
        
        % Get number trajectories already in group of current bin
        if isfield(bins(j), 'group')
            num_schedules = size(bins(j).group,1);  % current number of schedules in that bin
        else
            num_schedules = 0;
            bins(j).group = {};
        end
        
        % Copy over the trajectory
        bins(j).group(num_schedules+1,:) = trajectories(i,:);  % copy in that row to cell array
    
    end

end