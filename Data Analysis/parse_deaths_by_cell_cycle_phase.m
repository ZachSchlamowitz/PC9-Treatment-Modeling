%% parse_deaths_by_cell_cycle_phase.m
% This script determines what fraction of deaths at each timepoint occur in
% cells in G1 phase versus S/G2 phases of the cell cycle, to assess whether
% G1 cells are more likely to die than S/G2 cells (for pc9 cells).

% Author: Zach Schlamowitz (8/13/21)

%% Main Algorithm
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
        if 
        
    