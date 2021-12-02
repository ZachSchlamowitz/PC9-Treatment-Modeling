%% HML_dose_simulation_1drug.m
% Script to simulate different treatment plans for a 1-drug, 3-day,
% daily-dose treatment of PC9s with an EGFR inhibitor (e.g., osimertinib).
% Given set values of a "high" dose (H), "medium" dose (M), and "low" dose 
% (L), the script simulates cell growth according to that treatment plan
% using Pc9_ODE_Model_1drugsim.m, the 1 Drug ODE Model of PC-9 Cell 
% Response to treatment with varying doses of EGFR inhibition. The script 
% outputs a mapping of each (biologically) feasible treatmant schedule
% to the remaining number of cells, scaled by the total amount of drug
% adminstered.

% Author: Zach Schlamowitz (8/11/21)

%% CODE
% Initialize doses (nM)
H = 1.150;%1250;
M = 0.650;%600;
L = 0.200;%250;

% Initialize MTD
MTD = 1.45669; % µM

% Test to make sure dose choices fit all assumptions/conditions
assert(H*(.5^.5) + H > MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(M*(.5^.5) + H > MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(H*(.5^.5) + M > MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(H*(.5^.5) + L < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(L*(.5^.5) + H < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert((M*(.5^.5) + M)*(.5^.5) + M < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert((H*(.5^.5) + L)*(.5^.5) + M < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")

% Set simulation length
num_days = 7;

% Get all possible permutations of doses (i.e. theoretical schedules)
tic
fprintf('Obtaining all possible dosing schedules. \n')
dose_schedules = get_dosing_schedules(num_days, H, M, L);
toc

% Get rid of biologically infeasible schedules (those that exceed MTD)

for i = 1:size(dose_schedules,1)
    remainder = dose_schedules(i,1);
    for j = 1:size(dose_schedules,2)-1
        % Find dose after current day and add next day's treatment
        remainder = get_24hr_remaining_dose(remainder); % NOTE: Assume drug half life of 48hrs here, for consistency
        remainder = remainder + dose_schedules(i,j+1);

        % Check if feasible; replace that schedule with NaNs if not
        if remainder > MTD 
            dose_schedules(i,:) = nan(1,size(dose_schedules,2));
            break
        end
    end
end

%% Simulate ODE model for each valid combination of doses
% Max population over simulation period (i.e., total population when dose strategy is 0µM every day)
max_pop = get_max_pop(num_days); % obtained by running simulation with all zeros as [drug] inputs (for 3 days: 747.461853668818)

trajectories = cell(1+size(dose_schedules,1),7);  % initialize output cell array

for i = 1:size(dose_schedules,1)  % loop over dose schedules
    if ~isnan(dose_schedules(i,:)) % only consider valid dose schedules
        [pop_trajectory, drug_trajectory] = variable_daily_dose_sim(dose_schedules(i,:));  % simulate dose schedule's effects on cell populations
        
        % Store outputs
        trajectories{i,1} = get_dose_name(dose_schedules(i,:), H, M, L);
        trajectories(i,2:3) = {pop_trajectory, drug_trajectory};
        trajectories{i,4} = pop_trajectory(end,3);  % total population remaining at end of simulation
        trajectories{i,5} = max_pop - pop_trajectory(end,3);  % total cell death achieved by dose schedule
        trajectories{i,6} = sum(dose_schedules(i,:));  % Total drug administered
        trajectories{i,7} = trajectories{i,5}/trajectories{i,6};  % Effective Death: total cell death achieved by dose schedule normalized by total amount of drug given (i.e., number of dead cells per unit drug, µM)
    end
end

% To identify which dosing strategy was most effective, we rank strategies by their Death Counts:
% First we must remove empty row entries, storing good rows into new cell array "cleaned_trajectories"
cleaned_trajectories = {};
j=1;
for i=1:size(trajectories,1)
    if ~isempty(trajectories{i,1})
        cleaned_trajectories(j,:) = trajectories(i,:);
        j = j+1;
    end
end

% Sort by death count
[ranked_effdeaths, indices] = sort(cell2mat(cleaned_trajectories(:,5)));
sorted_trajectories = cleaned_trajectories(indices,:);

% Add Titles
% (Make space first)
trajectories(2:1+size(trajectories,1),:) = trajectories;
sorted_trajectories(2:1+size(sorted_trajectories,1),:) = sorted_trajectories;

trajectories{1,1} = "Dose Schedule";
trajectories{1,2} = "Pop Trajectories";
trajectories{1,3} = "Drug Trajectory";
trajectories{1,4} = "Final Total Pop";
trajectories{1,5} = "Death Count";
trajectories{1,6} = "Total Drug";
trajectories{1,7} = "Death /unit Drug";
sorted_trajectories{1,1} = "Dose Schedule";
sorted_trajectories{1,2} = "Pop Trajectories";
sorted_trajectories{1,3} = "Drug Trajectory";
sorted_trajectories{1,4} = "Final Total Pop";
sorted_trajectories{1,5} = "Death Count";
sorted_trajectories{1,6} = "Total Drug";
sorted_trajectories{1,7} = "Death /unit Drug";

%% Aux Functions
function remainder = get_24hr_remaining_dose(initial_dose)
    drug_half_life = 48; % (hrs) default: 48
    remainder = initial_dose * (1/2)^(24/drug_half_life);  % FLAG half life
end

function name = get_dose_name(dose_strategy, H, M, L)
    name = "";
    for j = 1:size(dose_strategy,2)
        if dose_strategy(j) == H
            name = name + "H";
        elseif dose_strategy(j) == M
            name = name + "M";
        elseif dose_strategy(j) == L
            name = name + "L";
        end
    end      
end

function max_pop = get_max_pop(num_days)
    [pop_trajectory, drug_trajectory] = variable_daily_dose_sim(zeros(1,num_days));
    max_pop = pop_trajectory(end,3);
end
