%% get_death_div_stats.m
%
% NOTE: This is a first pass attempt (succesful) at obtaining these stats,
% and uses an "average" algorithm approach that looks at the entire 
% population of cells over the 72hr period to pull out per hour values. See
% get_div_time_gefitnib.m for a more sophisticated approach. (8/2/21)
%
% This script obtains statistics regarding the death and division 
% rates/probabilities of cells (e.g., pc-9 non-small cell lung cancer
% cells) by analyzing their single-cell trace data. This data is to be 
% inputted in a matlab struct via a command:
% e.g.) >>    load('CleanGef0_1_200cells');
% Furthermore, it should then be added as a column to the struct matrix
% 'datasets':
% e.g.) >>    datasets = [MyOtherData CleanGef0_1_200cells Etc.];
% Note that the name of the struct variable that goes in this matrix may
% not default to the name of the .mat file which houses it. Verify this by
% opening the file manually. 

% Finally, it should also be noted that this script outputs statistics for
% both 72 and 24 hour segments of a 72 hour movie. However, to ensure
% accurate data, we ignore any cells with an odd number of divisions in the
% time period. This means that an overly high number of cells are ignored
% in the 24 hour case, as many have begun divisions but not completed them
% at hour 24. This may skew that statistic somewhat; keep this in mind.

% To only run death or division stats, comment out the other, undesired for
% loop. Or, if you're really feeling it, add in an if-else settings switch.

% Author: Zach Schlamowitz, 7/27/21

%% Initialization
% Initialize Output Structs
gef_stats_0p1uM = struct();
gef_stats_0p5uM = struct();
gef_stats_1uM = struct();

% Get dataset structs and arrange in one struct
load('CleanGef0_1_200cells');
load('Clean_Gef0_05_174cells');
load('Clean_Gef_1uM_200cells');

datasets = [CleanGef0_1 Clean_Gef0_05_200cells Clean_Gef_1uM_200cells]; % [CleanGef0_1] [Clean_Gef0_05_200cells] [Clean_Gef_1uM_200cells]


%% Death
% Loop over relevant datasets and run statistic extraction algorithm
for k=1:numel(datasets)
    
    full_72 = struct('deathMatrix',datasets(k).deathMatrixDataset);
    first_24 = struct('deathMatrix',datasets(k).deathMatrixDataset(:,1:(24*4)));  % Death in first 24 hours
    data_subsets = [full_72 first_24];
    subset_stats = [struct() struct()];
    
    for j=1:numel(data_subsets)
        num_deaths_per_cell = sum(data_subsets(j).deathMatrix,2);
        death_flags = zeros(numel(num_deaths_per_cell),1);
        invalid_deaths = 0;
        for i=1:numel(num_deaths_per_cell)
            if num_deaths_per_cell(i) > 1
                fprintf('death flag cell ')
                fprintf(num2str(i))
                fprintf('\n')
                death_flags(i) = 1;
                invalid_deaths = invalid_deaths + num_deaths_per_cell(i);
            end
        end
        % Store stats computed for this sub-dataset (either first 24 or all
        % 72 hrs) of this dose of gefitnib
        subset_stats(j).num_deaths_per_cell = num_deaths_per_cell;
        subset_stats(j).death_flags = death_flags;
        subset_stats(j).invalid_deaths = invalid_deaths;
    
    % Visual spacer
    if invalid_deaths > 0
        fprintf('\n')   
    end
    
    end
    
    total_deaths_72hrs = sum(subset_stats(1).num_deaths_per_cell) - subset_stats(1).invalid_deaths;
    valid_death_cells_72hrs = numel(subset_stats(1).num_deaths_per_cell) - sum(subset_stats(1).death_flags);
    prob_death_72hrs = total_deaths_72hrs / valid_death_cells_72hrs;
    
    total_divs_24hrs = sum(subset_stats(2).num_deaths_per_cell) - subset_stats(2).invalid_deaths;
    valid_death_cells_24hrs = numel(subset_stats(2).num_deaths_per_cell) - sum(subset_stats(2).death_flags);
    prob_death_24hrs = total_divs_24hrs / valid_death_cells_24hrs;
    
    output_structs(k).num_cells = numel(subset_stats(1).num_deaths_per_cell);
    output_structs(k).death_events_per_cell = subset_stats(1).num_deaths_per_cell;
    output_structs(k).death_flags_72hrs = subset_stats(1).death_flags;
    output_structs(k).num_invalid_deaths_72hrs = subset_stats(1).invalid_deaths;
    output_structs(k).num_cells_good4death_72hrs = valid_death_cells_72hrs;
    output_structs(k).num_deaths_72hrs = total_deaths_72hrs;
    output_structs(k).prob_death_72hrs =  prob_death_72hrs;
    output_structs(k).death_flags_24hrs = subset_stats(2).death_flags;
    output_structs(k).num_invalid_deaths_24hrs = subset_stats(2).invalid_deaths;
    output_structs(k).num_cells_good4death_24hrs = valid_death_cells_24hrs;
    output_structs(k).num_deaths_24hrs = total_divs_24hrs;
    output_structs(k).prob_death_24hrs =  prob_death_24hrs;
end


%% Division
for k=1:numel(datasets)
    
    full_72 = struct('divisionMatrix',datasets(k).divisionMatrixDataset);
    first_24 = struct('divisionMatrix',datasets(k).divisionMatrixDataset(:,1:(24*4)));  % Death in first 24 hours
    data_subsets = [full_72 first_24];
    subset_stats = [struct() struct()];
    
    for j=1:numel(data_subsets)
        num_div_events_per_cell = sum(data_subsets(j).divisionMatrix,2);
        div_flags = zeros(numel(num_div_events_per_cell),1);
        invalid_div_events = 0;
        for i=1:numel(num_div_events_per_cell)
            if mod(num_div_events_per_cell(i),2) ~= 0 % odd number of div events
                fprintf('div flag cell ')
                fprintf(num2str(i))
                fprintf('\n')
                div_flags(i) = 1;
                invalid_div_events = invalid_div_events + num_div_events_per_cell(i);
            end
        end
        % Store stats computed for this sub-dataset (either first 24 or all
        % 72 hrs) of this dose of gefitnib
        subset_stats(j).num_div_events_per_cell = num_div_events_per_cell;
        subset_stats(j).div_flags = div_flags;
        subset_stats(j).invalid_div_events = invalid_div_events;
     
    % Visual spacer
    if invalid_div_events > 0
        fprintf('\n')   
    end
    
    end
    
    total_divs_72hrs = (sum(subset_stats(1).num_div_events_per_cell) - subset_stats(1).invalid_div_events) / 2; % scalar
    valid_div_cells_72hrs = numel(subset_stats(1).num_div_events_per_cell) - sum(subset_stats(1).div_flags); % scalar
    avg_divs_per_cell_72 = total_divs_72hrs / valid_div_cells_72hrs; % scalar
    avg_divs_per_cell_per_hr_72 = avg_divs_per_cell_72/72; % scalar
    hrs_per_div_72 = 1/avg_divs_per_cell_per_hr_72; % scalar; number of hours required for one cell to copmlete one division (one cell cycle)
    
    total_divs_24hrs = (sum(subset_stats(2).num_div_events_per_cell) - subset_stats(2).invalid_div_events) / 2;
    valid_div_cells_24hrs = numel(subset_stats(2).num_div_events_per_cell) - sum(subset_stats(2).div_flags);
    avg_divs_per_cell_24 = total_divs_24hrs / valid_div_cells_24hrs; % scalar
    avg_divs_per_cell_per_hr_24 = avg_divs_per_cell_24/24; % scalar
    hrs_per_div_24 = 1/avg_divs_per_cell_per_hr_24; % scalar; number of hours required for one cell to copmlete one division (one cell cycle) in first 24 hours
    
    output_structs(k).num_cells = numel(subset_stats(1).num_div_events_per_cell);
    output_structs(k).div_events_per_cell = subset_stats(1).num_div_events_per_cell;
    output_structs(k).div_flags_72hrs = subset_stats(1).div_flags;
    output_structs(k).num_invalid_div_events_72hrs = subset_stats(1).invalid_div_events;
    output_structs(k).num_cells_good4divs_72hrs = valid_div_cells_72hrs;
    output_structs(k).num_divs_72hrs = total_divs_72hrs;
    output_structs(k).hrs_per_div_72 = hrs_per_div_72;
    output_structs(k).div_flags_24hrs = subset_stats(2).div_flags;
    output_structs(k).num_invalid_div_events_24hrs = subset_stats(2).invalid_div_events;
    output_structs(k).num_cells_good4divs_24hrs = valid_div_cells_24hrs;
    output_structs(k).num_divs_24hrs = total_divs_24hrs;
    output_structs(k).hrs_per_div_24 = hrs_per_div_24;

end