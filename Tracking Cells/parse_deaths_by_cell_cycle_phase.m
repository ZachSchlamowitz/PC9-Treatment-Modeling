%% parse_deaths_by_cell_cycle_phase.m
% This script determines what fraction of deaths at each timepoint occur in
% cells in G1 phase versus S/G2 phases of the cell cycle, to assess whether
% G1 cells are more likely to die than S/G2 cells (for pc9 cells).

% Author: Zach Schlamowitz (8/13/21)

%% Clean Up Data
load('fucci_gefitnib_1uM_RFP');
load('fucci_gefitnib_1uM_GFP');
RFP_dataset = fucci_gefitnib_1uM_RFP;
GFP_dataset = fucci_gefitnib_1uM_GFP;

[death_RFP, div_RFP, intensity_RFP] = cleanup_data(RFP_dataset);
[death_GFP, div_GFP, intensity_GFP] = cleanup_data(GFP_dataset);

%% Main Algorithm
num_cells = size(death_RFP,1); % number of cells tracks in data file
num_timepoints = 282;

% Channel intensity thresholds to determine if a cell is in G1 or S/G2.
% These are obtained by inspection (looking at images to approximate a
% threshold value).
G1_threshold = 2400;
S_G2_threshold = 3100;

%%% Initialize output matrices %%%
% Cell Counts
num_G1_cells = zeros(1,282);
num_SG2_cells = zeros(1,282);
% num_bothon_cells = zeros(1,282);
% num_bothoff_cells = zeros(1,282);
% num_nan = zeros(1,282);
num_other_cells = zeros(1,282);

% Death Frac Type I: fraction of each population that dies
frac_G1_cells_dying = zeros(1,282);
frac_SG2_cells_dying = zeros(1,282);
frac_other_cells_dying = zeros(1,282);

% Cumulative Value of Death Frac Type I
cumfrac_G1_cells_dying = zeros(1,282);
cumfrac_SG2_cells_dying = zeros(1,282);
cumfrac_other_cells_dying = zeros(1,282);

% Death Counts
G1_deaths = zeros(1,282);
G1_deaths_cumulative = zeros(1,282);
SG2_deaths = zeros(1,282);
SG2_deaths_cumulative = zeros(1,282);
other_deaths = zeros(1,282);
other_deaths_cumulative = zeros(1,282);

%%% Loop over each timepoint for every cell and categorize it by cell cycle phase;
% categorize deaths by doing same for death matrix %%%
for i = 1:num_cells
    for j = 1:num_timepoints
        
        % Categorize the cell (at the current timepoint)
        if intensity_RFP(i,j) >= G1_threshold && intensity_GFP(i,j) < S_G2_threshold
            num_G1_cells(j) = num_G1_cells(j) + 1;                
        elseif intensity_RFP(i,j) < G1_threshold && intensity_GFP(i,j) >= S_G2_threshold
            num_SG2_cells(j) = num_SG2_cells(j) + 1;
%         elseif intensity_RFP(i,j) < G1_threshold && intensity_GFP(i,j) < S_G2_threshold
%             num_bothoff_cells(j) = num_bothoff_cells(j) + 1;
%         elseif intensity_RFP(i,j) >= G1_threshold && intensity_GFP(i,j) >= S_G2_threshold
%             num_bothon_cells(j) = num_bothon_cells(j) + 1;
        elseif isnan(intensity_RFP(i,j)) || isnan(intensity_GFP(i,j)) || intensity_RFP(i,j) == -1 || intensity_GFP(i,j) == -1
            fprintf('') % do nothing
        else
            num_other_cells(j) = num_other_cells(j) + 1;
        end
        
        % Categorize deaths (at current timepoint by looking at signal of previous timepoint)
        if death_RFP(i,j) == 1                   
            if intensity_RFP(i,j-1) >= G1_threshold && intensity_GFP(i,j-1) < S_G2_threshold
                G1_deaths(j) = G1_deaths(j) + 1;
            elseif intensity_RFP(i,j-1) < G1_threshold && intensity_GFP(i,j-1) >= S_G2_threshold
                SG2_deaths(j) = SG2_deaths(j) + 1;
            elseif isnan(intensity_RFP(i,j-1)) || isnan(intensity_GFP(i,j-1)) || intensity_RFP(i,j-1) == -1 || intensity_GFP(i,j-1) == -1
                fprintf('Problematic categorization. Skipping this frame for this cell.\n')
                continue
            else
                other_deaths(j) = other_deaths(j) + 1;
            end            
        end
                
        % Record death counts and fracs
        frac_G1_cells_dying(j) = G1_deaths(j)/num_G1_cells(j);
        frac_SG2_cells_dying(j) = SG2_deaths(j)/num_SG2_cells(j);
        frac_other_cells_dying(j) = other_deaths(j)/num_other_cells(j);

        cumfrac_G1_cells_dying(j) = sum(frac_G1_cells_dying(1:j));
        cumfrac_SG2_cells_dying(j) = sum(frac_SG2_cells_dying(1:j));
        cumfrac_other_cells_dying(j) = sum(frac_other_cells_dying(1:j));
        
        G1_deaths_cumulative(j) = sum(G1_deaths(1:j));
        SG2_deaths_cumulative(j) = sum(SG2_deaths(1:j));
        other_deaths_cumulative(j) = sum(other_deaths(1:j));
        
    end
end

total_deaths_cumulative = G1_deaths_cumulative + SG2_deaths_cumulative + other_deaths_cumulative;

% Death Frac Type II: fraction of deaths which occur in each cell cycle phase
cumfrac_deaths_in_G1 = zeros(1,282);
cumfrac_deaths_in_SG2 = zeros(1,282);
cumfrac_deaths_in_other = zeros(1,282);

for j = 1:size(total_deaths_cumulative,2)
    if total_deaths_cumulative(j) == 0
        continue
    else
        cumfrac_deaths_in_G1(j) = G1_deaths_cumulative(j)/total_deaths_cumulative(j);
        cumfrac_deaths_in_SG2(j) = SG2_deaths_cumulative(j)/total_deaths_cumulative(j);
        cumfrac_deaths_in_other(j) = other_deaths_cumulative(j)/total_deaths_cumulative(j);

    end
end

categorized_num_cells = [num_G1_cells; num_SG2_cells; num_other_cells]; %num_bothon_cells; num_bothoff_cells];%; num_other_cells; num_nan];

categorized_deaths = [G1_deaths; SG2_deaths; other_deaths];
categorized_deaths_cumulative = [G1_deaths_cumulative; SG2_deaths_cumulative; other_deaths_cumulative];
categorized_death_fracs_T1 = [frac_G1_cells_dying; frac_SG2_cells_dying; frac_other_cells_dying];
categorized_death_fracs_T1_cumulative = [cumfrac_G1_cells_dying; cumfrac_SG2_cells_dying; cumfrac_other_cells_dying];
categorized_death_fracs_T2_cumulative = [cumfrac_deaths_in_G1; cumfrac_deaths_in_SG2; cumfrac_deaths_in_other];

%% %% %% Plotting %% %% %%
%%% Calls %%%
plot_num_cells(categorized_num_cells);
plot_death_frac_T1(categorized_death_fracs_T1);
plot_death_frac_T1_cum(categorized_death_fracs_T1_cumulative);
plot_death_frac_T2_cum(categorized_death_fracs_T2_cumulative);
plot_cum_deaths(categorized_deaths_cumulative);
% static_death_fracs_T2_cum(categorized_death_fracs_0p1uM, categorized_death_fracs_1uM)

%%% DEATH PLOT FUNCS %%%
function plot_death_frac_T1(categorized_death_fracs)
    f = figure;
    f.Position = [1100 400 1000 500];
    p = plot([1:282], categorized_death_fracs, '-');
    for i=1:3
        p(i).LineWidth = 3;
    end
    t = title('Fraction of G1 Cells and S/G2 Cells that Die per Timestep, (1uM Gefitinib)');
    t.FontSize = 20;

    x = xlabel('Time (hrs)');
    x.FontSize = 18;
    ticks = [0:8:282];
    labels = [0:2:72];
    xticks(ticks)
    xticklabels(labels)

    y = ylabel('Fraction of Subpopulation');
    y.FontSize = 18;

    l = legend('G1', 'S/G2', 'Failed to Distinguish', 'Location', 'northwest');
    l.FontSize = 15;
end
function plot_death_frac_T1_cum(categorized_death_fracs)
    f = figure;
    f.Position = [1100 400 1000 500];
    p = plot([1:282], categorized_death_fracs, '-');
    for i=1:3
        p(i).LineWidth = 3;
    end
    t = title('Cumulative Fraction of G1 Cells and S/G2 Cells that Die per Timestep (1uM Gefitinib)');
    t.FontSize = 20;

    x = xlabel('Time (hrs)');
    x.FontSize = 18;
    ticks = [0:8:282];
    labels = [0:2:72];
    xticks(ticks)
    xticklabels(labels)

    y = ylabel('Fraction of Subpopulation');
    y.FontSize = 18;

    l = legend('G1', 'S/G2', 'Failed to Distinguish', 'Location', 'northwest');
    l.FontSize = 15;
end

function plot_death_frac_T2_cum(categorized_death_fracs)
    f = figure;
    f.Position = [1100 400 1000 500];
    p = plot([1:282], categorized_death_fracs, '-');
    for i=1:3
        p(i).LineWidth = 3;
    end
    t = title('Cumulative Fraction of Deaths Which Occur in G1 vs S/G2 per Timestep, (1uM Gefitinib)');
    t.FontSize = 20;

    x = xlabel('Time (hrs)');
    x.FontSize = 18;
    ticks = [0:8:282];
    labels = [0:2:72];
    xticks(ticks)
    xticklabels(labels)

    y = ylabel('Fraction of All Deaths');
    y.FontSize = 18;

    l = legend('G1', 'S/G2', 'Failed to Distinguish', 'Location', 'northwest');
    l.FontSize = 15;
end

function plot_cum_deaths(categorized_deaths_cumulative)
    f = figure;
    f.Position = [50 400 1000 500];
    p = plot([1:282], categorized_deaths_cumulative, '-');
    for i=1:3
        p(i).LineWidth = 3;
    end
    t = title('Deaths per Timestep, Categorized by Cell Cycle Phase, (1uM Gefitinib)');
    t.FontSize = 20;

    x = xlabel('Time (hrs)');
    x.FontSize = 18;
    ticks = [0:8:282];
    labels = [0:2:72];
    xticks(ticks)
    xticklabels(labels)

    y = ylabel('Cumulative Deaths');
    y.FontSize = 18;

    l = legend('G1', 'S/G2', 'Failed to Distinguish', 'Location', 'northwest');
    l.FontSize = 15;
end

function static_death_fracs_T2_cum(categorized_death_fracs_0p1uM, categorized_death_fracs_1uM)
    % NOTE: need to have run the main algorithm twice and saved outputs to
    % input different dose data into this function!
    categorized_death_fracs_static = ...
       [categorized_death_fracs_0p1uM(1,282) categorized_death_fracs_1uM(1,282); 
        categorized_death_fracs_0p1uM(2,282) categorized_death_fracs_1uM(2,282); 
        categorized_death_fracs_0p1uM(3,282) categorized_death_fracs_1uM(3,282); 
        categorized_death_fracs_0p1uM(4,282) categorized_death_fracs_1uM(4,282)];
    
    f = figure;
    f.Position = [50 400 1000 500];
    xlabels = categorical({'G1' 'S' 'G2' 'Unknown'});
    b = bar(xlabels, categorized_death_fracs_static);
    t = title('Deaths Fractions after 72 hours (Categorized by Cell Cycle Phase and Gefitinib Dose)');
    t.FontSize = 20;
    ylim([0 1])

    x = xlabel('Cell Cycle Phase');
    x.FontSize = 18;

    y = ylabel('Fraction of Total Deaths');
    y.FontSize = 18;
    
    l = legend('0.1uM Gefitinib','1uM Gefitinib');
    l.FontSize = 15;
    
    % Display values above bars
    xtips1 = b(1).XEndPoints;
    ytips1 = b(1).YEndPoints;
    labels1 = string(b(1).YData);
    text(xtips1,ytips1,labels1,'HorizontalAlignment','center','VerticalAlignment','bottom')
    xtips2 = b(2).XEndPoints;
    ytips2 = b(2).YEndPoints;
    labels2 = string(b(2).YData);
    text(xtips2,ytips2,labels2,'HorizontalAlignment','center','VerticalAlignment','bottom')
    
end

%%% OTHER PLOT FUNCS %%%
function plot_num_cells(categorized_num_cells)
    f = figure;
    f.Position = [50 400 1000 500];
    p = plot([1:282], categorized_num_cells, '-');
    for i=1:size(p)
        p(i).LineWidth = 3;
    end
    t = title('Number of Cells in G1 vs S/G2 per Timestep (1uM Gefitinib)');
    t.FontSize = 20;

    x = xlabel('Time (hrs)');
    x.FontSize = 18;
    ticks = [0:8:282];
    labels = [0:2:72];
    xticks(ticks)
    xticklabels(labels)

    y = ylabel('Number of Cells');
    y.FontSize = 18;

%     l = legend('G1', 'S/G2', 'Both Reporters On', 'Both Reporters Off', 'Location', 'northwest');
    l = legend('G1', 'S/G2', 'Failed to Distinguish', 'Location', 'northwest');
    l.FontSize = 15;
end

%% Clean Up Function
function [death,div,intensity] = cleanup_data(dataset)
    % Extract death, division, and intensity matrices from dataset
    intensity_matrix = dataset(1).singleCellTraces;
    div_matrix = dataset(1).divisionMatrixDataset;
    death_matrix = dataset(1).deathMatrix;

    % Identify questionable cells in both division and death datasets
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
    intensity_matrix_cleaned = zeros(size(intensity_matrix));
    div_matrix_cleaned = zeros(size(div_matrix)); 
    death_matrix_cleaned = zeros(size(death_matrix));

    % Copy original data to cleaned matrices replacing questionable cells with rows of NaNs
    for i = 1:size(div_matrix,1)
        if death_flags(i) == 1 || div_flags(i) == 1
            intensity_matrix_cleaned(i,:) = NaN(1,size(intensity_matrix,2));
            death_matrix_cleaned(i,:) = NaN(1,size(death_matrix,2));
            div_matrix_cleaned(i,:) = NaN(1,size(div_matrix,2));
        else % No problem, just copy over old data
            intensity_matrix_cleaned(i,:) = intensity_matrix(i,:);
            death_matrix_cleaned(i,:) = death_matrix(i,:);        
            div_matrix_cleaned(i,:) = div_matrix(i,:);
        end
    end

    % Convert matrix entries for dead cells to NaNs
    % ex)    % death_matrix = [0 0 1 0;
             %                 0 1 0 0;
             %                 1 0 1 0];
             % out_data = zeros(3,4);

    for i = 1:size(death_matrix,1)
        for j = 1:size(death_matrix,2)
            if death_matrix(i,j) == 1
                for k = j:size(death_matrix,2) % go to the end of the row
                    div_matrix_cleaned(i,k) = NaN;
                    
                    try
                        assert((k+1)<=size(death_matrix,2))
                    catch
                        continue
                    end
                    death_matrix_cleaned(i,k+1) = NaN;
                    
                    try
                        assert(intensity_matrix_cleaned(i,k+1) == -1, 'Cell claimed to be dead but intensity entry not -1 at next entry')
                    catch
                        fprintf('i = ')
                        fprintf(num2str(i))
                        fprintf('  j = ')
                        fprintf(num2str(j))
                        fprintf('  next intensity value = ')
                        fprintf(num2str(intensity_matrix_cleaned(i,k+1)))
                        fprintf('\n')
                    end
                    intensity_matrix_cleaned(i,k) = NaN;
                end
                break  % advance to next row (i value)
            end
        end
    end

    death = death_matrix_cleaned;
    div = div_matrix_cleaned;
    intensity = intensity_matrix_cleaned;
    
    
end
