%% get_div_time_gefitnib.m
% 
% This script takes in a Division matrix from single cell trace data and
% returns the average percent of cells which divide in any given hour. This
% is determined by creating a matrix of the same dimensions in which dead
% cells are marked with NaN entries rather than zeros. Then, the first
% division marker in any pair of divisions is converted to a 0. The percent
% of cells with a 1 in any given column is then computed to provide the
% percent of the population which divides in each hour. Finally, this value
% is averaged. 

% Author: Zach Schlamowitz (7/29/21)

%% Initialization

% Run settings
ignore_first_division = 1;


% Get dataset struct and extract division matrix
% load('CleanGef0_1_200cells');
% dataset = CleanGef0_1;
% load('Clean_Gef0_05_174cells');
% dataset = Clean_Gef0_05_200cells;
load('Clean_Gef_1uM_200cells');
dataset = Clean_Gef_1uM_200cells;

in_data = dataset(1).divisionMatrixDataset;

death_matrix = dataset(1).deathMatrixDataset;

%% Remove Questionable cells from both division and death datasets
num_deaths_per_cell = sum(death_matrix,2);
death_flags = zeros(numel(num_deaths_per_cell),1);
invalid_deaths = 0;
num_div_events_per_cell = sum(in_data,2);
div_flags = zeros(numel(num_div_events_per_cell),1);
invalid_div_events = 0;
for i=1:numel(num_deaths_per_cell)
    if num_deaths_per_cell(i) > 1
        fprintf('death flag cell ')
        fprintf(num2str(i))
        fprintf('\n')
        death_flags(i) = 1;
        invalid_deaths = invalid_deaths + num_deaths_per_cell(i);
    end
    
    if mod(num_div_events_per_cell(i),2) ~= 0 % odd number of div events
        fprintf('div flag cell ')
        fprintf(num2str(i))
        fprintf('\n')
        div_flags(i) = 1;
        invalid_div_events = invalid_div_events + num_div_events_per_cell(i);
    end
end


out_data = []; % initialize matrix for analysis

% Copy in_data to out_data replacing questionable cells with rows of NaNs
for i = 1:size(in_data,1)
    if death_flags(i) == 1 || div_flags(i) == 1
        out_data(i,:) = NaN(1,size(in_data,2));
    else
        out_data(i,:) = in_data(i,:);
    end
end

%% Convert division matrix entries for dead cells to NaNs
% death_matrix = [0 0 1 0;
%                 0 1 0 0;
%                 1 0 1 0];
% out_data = zeros(3,4);
% 

for i = 1:size(death_matrix,1)
    for j = 1:size(death_matrix,2)
        if death_matrix(i,j) == 1
            for k = j:size(death_matrix,2) % go to the end of the row
                out_data(i,k) = NaN;
            end
            break  % advance to next row (i value)
        end
    end
end
      
%% Get rid of first 1 in pair of 1s indicating division
% out_data = [0 0 1 0 0 1 0;
%             0 1 0 1 0 1 1;
%             1 0 0 0 0 0 1];

for i = 1:size(out_data,1)
    div_count = 0;
    for j = 1:size(out_data,2)
        if out_data(i,j) == 1
            div_count = div_count + 1;
            if mod(div_count,2) == 1
                out_data(i,j) = 0;
            end
        end
    end
end
                
%% Compute percent of cells which divide in any hour
% out_data = [0   NaN 0 0 1 0 0;
%             NaN 0   0 1 0 0 1;
%             0   0   0 0 0 0 1];


num_divs_per_ts = zeros(1,size(out_data,2)); % NOTE: 'ts' is timestep
num_alive_cells_per_ts = zeros(1,size(out_data,2));
frac_divs_per_ts = zeros(1,size(out_data,2));

for j = 1:size(out_data,2)
    num_cells = 0;
    num_divs = 0;
    for i = 1:size(out_data,1)
        if isnan(out_data(i,j))
            continue
        end
        num_cells = num_cells + 1;
        num_divs = num_divs + out_data(i,j);
    end
    
    num_divs_per_ts(1,j) = num_divs;
    num_alive_cells_per_ts(1,j) = num_cells;
    frac_divs_per_ts(1,j) = num_divs/num_cells;
    
end

avg_frac_divs_in_1ts = mean(frac_divs_per_ts);
avg_frac_divs_in_1hr = avg_frac_divs_in_1ts * 4; % assumes 15min ts
hrs_to_div = 1/avg_frac_divs_in_1hr;

%% Compute underestimate value
% To account for the fact that many cells divide once, as normal, before 
% responding to the drug, we compute the same statistics, this time
% ignoring the first division for any cell. This provides an underestimate
% for the time to division, which, combined with the overestimate given by
% including the "drug-uninhibited" first division, provides bounds for the 
% real division time. 
% num_div_events_per_cell = [2; 0; 1;];
% out_data = [0 0 1 0 0 0 0;
%             0 0 0 0 0 0 0;
%             1 0 0 0 0 0 1];

% CONVERT ALL ENTRIES THROUGH THE FIRST DIVISION TO NaNs %
% if ignore_first_division == 1
%     for i = 1:size(out_data,1)
%         div_count = 0;
%         j = 1;
%         while div_count < 2            
%             if out_data(i,j) == 1
%                 div_count = div_count + 1;
%             end
%             
%             out_data(i,j) = NaN;
%             j = j + 1;
%         end
%     end     
% end


% CONVERT ALL ENTRIES THROUGH THE FIRST DIVISION TO NaNs %
if ignore_first_division == 1
    out_data_ifd = out_data;
    for i = 1:size(out_data_ifd,1)
        j = 1;
        
        if num_div_events_per_cell(i) == 0
            continue
        end
        
        has_divided = false;
        while ~has_divided          
            if out_data_ifd(i,j) == 1
                has_divided = true;
            end
            
            out_data_ifd(i,j) = NaN;
            
            if j == size(out_data_ifd,2)
                break
            else
                j = j + 1;
            end
            
        end
    end     
end

% Compute percent of cells which divide in any hour, now ignoring first
% divisions (ifd)
% out_data = [0   NaN 0 0 1 0 0;
%             NaN 0   0 1 0 0 1;
%             0   0   0 0 0 0 1];


num_divs_per_ts_ifd = zeros(1,size(out_data_ifd,2)); % NOTE: 'ts' is timestep
num_alive_cells_per_ts_ifd = zeros(1,size(out_data_ifd,2));
frac_divs_per_ts_ifd = zeros(1,size(out_data_ifd,2));

for j = 1:size(out_data_ifd,2)
    num_cells = 0;
    num_divs = 0;
    for i = 1:size(out_data_ifd,1)
        if isnan(out_data_ifd(i,j))
            continue
        end
        num_cells = num_cells + 1;
        num_divs = num_divs + out_data_ifd(i,j);
    end
    
    num_divs_per_ts_ifd(1,j) = num_divs;
    num_alive_cells_per_ts_ifd(1,j) = num_cells;
    frac_divs_per_ts_ifd(1,j) = num_divs/num_cells;
    
end

avg_frac_divs_in_1ts_ifd = mean(frac_divs_per_ts_ifd);
avg_frac_divs_in_1hr_ifd = avg_frac_divs_in_1ts_ifd * 4; % assumes 15min ts
hrs_to_div_ifd = 1/avg_frac_divs_in_1hr_ifd;

%% Convert hrs to div to hrs in G1 vs S/G2/M and compute transition rates for cell cycle
time_in_SG2M = 13;  % assumes 13 (hrs) spent in S/G2/M regardless of drug (based on generic 24hr cell cycle)
time_in_G1 = hrs_to_div - time_in_SG2M; % hrs
rate_of_S_entry = 1/time_in_G1; % fraction of G1 population which exits and enters S/G2/M per hour (alpha in model)
rate_of_M_entry = 1/time_in_SG2M; % fraction of S/G2/M cells which exit per hour (beta in model)

% To add these to new_gef_stats_____uM.mat files, use:
%     hrs_to_div = new_gef_stats_1uM.hrs_to_div;
% 
%     time_in_SG2M = 13;  % assumes 13 (hrs) spent in S/G2/M regardless of drug (based on generic 24hr cell cycle)
%     time_in_G1 = hrs_to_div - time_in_SG2M; % hrs
%     rate_of_S_entry = 1/time_in_G1; % fraction of G1 population which exits and enters S/G2/M per hour
%     rate_of_M_entry = 1/time_in_SG2M; % fraction of S/G2/M cells which exit per hour
% 
%     new_gef_stats_1uM.rate_of_S_entry = rate_of_S_entry;
%     new_gef_stats_1uM.rate_of_M_entry = rate_of_M_entry;
% 
%     clear hrs_to_div time_in_G1 rate_of_S_entry rate_of_M_entry