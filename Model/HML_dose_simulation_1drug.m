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
H = 1150;%1250;
M = 650;%600;
L = 200;%250;

% Initialize MTD
MTD = 1456.69; % nM

% Test to make sure dose choices fit all assumptions/conditions
assert(H*(.5^.5) + H > MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(M*(.5^.5) + H > MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(H*(.5^.5) + M > MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(H*(.5^.5) + L < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert(L*(.5^.5) + H < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert((M*(.5^.5) + M)*(.5^.5) + M < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")
assert((H*(.5^.5) + L)*(.5^.5) + M < MTD, "Chosen dose values (High, Medium, and Low) fail required conditions.")

% Get all possible combinations of doses (i.e. theoretical schedules)
dose_schedules = zeros(3^3, 3);
dose_schedules = [H H H; 
                  H H M; 
                  H H L; 
                  H M H; 
                  H M M; 
                  H M L;
                  H L H; 
                  H L M; %
                  H L L; %
                  
                  M H H; 
                  M H M; 
                  M H L; 
                  M M H; 
                  M M M; %
                  M M L; %
                  M L H; % FLAG: At H=1250,M=600,L=250 this does not work (=1727). But should it?
                  M L M; %
                  M L L; %
                  
                  L H H; 
                  L H M; 
                  L H L; %
                  L M H; 
                  L M M; %
                  L M L; %
                  L L H; % FLAG: At H=1250,M=600,L=250 this does not work (=1552) but at 1150/650/200 it does.
                  L L M; %
                  L L L];% 
              
 % Get rid of biologically infeasible schedules (those that exceed MTD)
for i = 1:size(dose_schedules,1)
   % Dose after day 2 and after day 3
   day2_remainder = dose_schedules(i,1)*(1/2)^(24/48) + dose_schedules(i,2);
   day3_remainder =    day2_remainder  *(1/2)^(24/48) + dose_schedules(i,3);
   if day2_remainder > MTD || day3_remainder > MTD
       dose_schedules(i,1:3) = [NaN NaN NaN];
   end
    
end

% Simulate ODE model for each valid combination of doses
trajectories = cell(size(dose_schedules,1),2);

for i = 1:size(dose_schedules,1)
    if dose_schedules(i,:) ~= [NaN NaN NaN]
        [pop_trajectory drug_trajectory] = variable_daily_dose_sim(dose_schedules(i,:));
        trajectories(i,:) = {pop_trajectory, drug_trajectory};
    end
end
