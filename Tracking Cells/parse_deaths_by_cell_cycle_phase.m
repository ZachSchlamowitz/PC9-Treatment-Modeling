%% parse_deaths_by_cell_cycle_phase.m
% This script determines what fraction of deaths at each timepoint occur in
% cells in G1 phase versus S/G2 phases of the cell cycle, to assess whether
% G1 cells are more likely to die than S/G2 cells (for pc9 cells).

% Author: Zach Schlamowitz (8/13/21)

%% Clean Up Data
% Get dataset struct and extract division matrix
load('CleanGef0_1_200cells');
dataset = CleanGef0_1;

intensity_matrix = dataset(1).singleCellTraces;
div_matrix = dataset(1).divisionMatrixDataset;
death_matrix = dataset(1).deathMatrixDataset;

% Remove Questionable cells from both division and death datasets
num_deaths_per_cell = sum(death_matrix,2);
death_flags = zeros(numel(num_deaths_per_cell),1);
invalid_deaths = 0;
num_div_events_per_cell = sum(div_matrix,2);
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

% Initialize matrices for analysis
intensity_matrix_cleaned = [];
div_matrix_cleaned = []; 
death_matrix_cleaned = [];

% Copy divisions to div_matrix_cleaned replacing questionable cells with rows of NaNs
for i = 1:size(div_matrix,1)
    if death_flags(i) == 1 || div_flags(i) == 1
        intensity_matrix_cleaned(i,:) = NaN(1,size(intensity_matrix,2));
        death_matrix_cleaned(i,:) = NaN((1,size(death_matrix,2));
        div_matrix_cleaned(i,:) = NaN(1,size(div_matrix,2));
    else
        intensity_matrix_cleaned(i,:) = intensity_matrix(i,:);
        death_matrix_cleaned(i,:) = death_matrix(i,:);        
        div_matrix_cleaned(i,:) = div_matrix(i,:);
    end
end

% Convert division matrix entries for dead cells to NaNs
    % death_matrix = [0 0 1 0;
    %                 0 1 0 0;
    %                 1 0 1 0];
    % out_data = zeros(3,4);

for i = 1:size(death_matrix,1)
    for j = 1:size(death_matrix,2)
        if death_matrix(i,j) == 1
            for k = j:size(death_matrix,2) % go to the end of the row
                div_matrix_cleaned(i,k) = NaN;
                death_matrix_cleaned(i,k) = NaN;
            end
            break  % advance to next row (i value)
        end
    end
end



pause

%% Main Algorithm
load('G1_intensity_matrix');
load('S_G2_intensity_matrix');
load('death_matrix')

num_cells = 51; % number of cells tracks in data file
num_timepoints = 282;

% Channel intensity thresholds to determine if a cell is in G1 or S/G2.
% These are obtained by inspection (looking at images to approximate a
% threshold value).
G1_threshold = 0;
S_G2_threshold = 0;

% Initialize output matrices
G1_deaths = zeros(1,282);
S_deaths = zeros(1,282);
G2_deaths = zeros(1,282);
other_deaths = zeros(1,282);

% Loop over entries of death matrix looking for death events; for each, 
% categorize it by cell cycle phase
for i = 1:num_cells
    for j = 1:num_timepoints
        if death_matrix(i,j) == 1
            if G1_intensity_matrix(i,j) >= G1_threshold && S_G2_intensity_matrix < S_G2_threshold
                G1_deaths(j) = G1_deaths(j) + 1;
            elseif G1_intensity_matrix(i,j) >= G1_threshold && S_G2_intensity_matrix >= S_G2_threshold
                S_deaths(j) = S_deaths(j) + 1;
            elseif G1_intensity_matrix(i,j) < G1_threshold && S_G2_intensity_matrix >= S_G2_threshold
                G2_deaths(j) = G2_deaths(j) + 1;
            else
                other_deaths(j) = other_deaths(j) + 1;
            end
        end
    end
end
    