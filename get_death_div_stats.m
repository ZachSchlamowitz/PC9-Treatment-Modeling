
% Initialize Output Structs
% gef_stats_0p1uM = struct();
gef_stats_0p5uM = struct();
% gef_stats_1uM = struct();

% Get dataset structs and arrange in one struct
% load('CleanGef0_1_200cells');
load('Clean_Gef0_05_174cells');
% load('Clean_Gef_1uM_200cells');

datasets = [Clean_Gef0_05_200cells]; % [CleanGef0_1]

% % For 1uM Gefitnib:
% num_divs_per_cell = sum(Clean_Gef_1uM_200cells.divisionMatrixDataset,2);
% div_flags = zeros(numel(num_divs_per_cell),1);
% for i=1:numel(num_divs_per_cell)
%     if mod(num_divs_per_cell(i),2) ~= 0
%         fprintf('flag cell ')
%         fprintf(num2str(i))
%         fprintf('\n')
%         div_flags(i) = 1;
%     end
% end
% 
% total_divs = sum(num_divs_per_cell);
% for i=[1 18 22 26 41 48 66 125 126]
%     total_divs = total_divs - num_divs_per_cell(i);
% end
% 
% gef_stats_1uM.div_events_per_cell = num_divs_per_cell;
% gef_stats_1uM.div_flags = div_flags;
% gef_stats_1uM.num_divs = total_divs/2;
% gef_stats_1uM.prob_div = gef_stats_1uM.num_divs/gef_stats_1uM.num_cells_good4divs
% 
% % For 0.1uM Gefitnib:
% num_divs_per_cell = sum(CleanGef0_1.divisionMatrixDataset,2);
% div_flags = zeros(numel(num_divs_per_cell),1);
% for i=1:numel(num_divs_per_cell)
%     if mod(num_divs_per_cell(i),2) ~= 0
%         fprintf('flag cell ')
%         fprintf(num2str(i))
%         fprintf('\n')
%         div_flags(i) = 1;
%     end
% end
% 
% total_divs = sum(num_divs_per_cell);
% for i=[33 67 68 134 165]
%     total_divs = total_divs - num_divs_per_cell(i);
% end
% 
% gef_stats_0p1uM.div_events_per_cell = num_divs_per_cell;
% gef_stats_0p1uM.div_flags = div_flags;
% gef_stats_0p1uM.num_divs = total_divs/2;
% gef_stats_0p1uM.prob_div = gef_stats_0p1uM.num_divs/gef_stats_0p1uM.num_cells_good4divs
% 
% For 0.5uM Gefitnib:
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
                fprintf('flag cell ')
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
        
    end

    total_deaths_72hrs = sum(subset_stats(1).num_deaths_per_cell) - subset_stats(1).invalid_deaths;
    valid_death_cells_72hrs = numel(subset_stats(1).num_deaths_per_cell) - sum(subset_stats(1).death_flags);
    prob_death_72hrs = total_deaths_72hrs / valid_death_cells_72hrs;
    
    total_deaths_24hrs = sum(subset_stats(2).num_deaths_per_cell) - subset_stats(2).invalid_deaths;
    valid_death_cells_24hrs = numel(subset_stats(2).num_deaths_per_cell) - sum(subset_stats(2).death_flags);
    prob_death_24hrs = total_deaths_24hrs / valid_death_cells_24hrs;
    
    output_structs(k).num_cells = numel(subset_stats(1).num_deaths_per_cell);
    output_structs(k).death_events_per_cell = subset_stats(1).num_deaths_per_cell;
    output_structs(k).num_cells_good4death_72hrs = valid_death_cells_72hrs;
    output_structs(k).num_deaths_72hrs = total_deaths_72hrs;
    output_structs(k).prob_death_72hrs =  prob_death_72hrs;
    output_structs(k).num_cells_good4death_24hrs = valid_death_cells_24hrs;
    output_structs(k).num_deaths_24hrs = total_deaths_24hrs;
    output_structs(k).prob_death_24hrs =  prob_death_24hrs;
end


% %% Division
% num_divs_per_cell = sum(Clean_Gef0_05_200cells.divisionMatrixDataset,2);
% div_flags = zeros(numel(num_divs_per_cell),1);
% for i=1:numel(num_divs_per_cell)
%     if mod(num_divs_per_cell(i),2) ~= 0
%         fprintf('flag cell ')
%         fprintf(num2str(i))
%         fprintf('\n')
%         div_flags(i) = 1;
%     end
% end
% 
% total_divs = sum(num_divs_per_cell);
% for i=[33 67 68 134 165]
%     total_divs = total_divs - num_divs_per_cell(i);
% end
% 
% gef_stats_0p1uM.div_events_per_cell = num_divs_per_cell;
% gef_stats_0p1uM.div_flags = div_flags;
% gef_stats_0p1uM.num_divs = total_divs/2;
% gef_stats_0p1uM.prob_div = gef_stats_0p1uM.num_divs/gef_stats_0p1uM.num_cells_good4divs
% 



