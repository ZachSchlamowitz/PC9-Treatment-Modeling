function [population_trajectories, drug_trajectory] = variable_daily_dose_sim(dose_by_day)

%% Pc9_ODE_Model_1drugsim.m
% 1 Drug ODE Model of PC-9 Cell Response to treatment with varying doses of EGFR inhibition
% Author: Zach Schlamowitz (8/11/21)

%% Main equations

% User Initializations:
num_days = size(dose_by_day, 2);
initconds = [40; 60; 0];  % number of cells in [G1 S/G2 XX] and EGFRi concentration [XX XX EGFRi]

% Internal Initializations:
% tspan = [0 72];
tspan = [0:24];  % Duration of one "simulation period" (time between dosing) in hours (default: [0:24] specifies 1 day)
population_trajectories = zeros(num_days*24,3);  % vector to house full trajectory of Vulnerable, Resistant, and Total populations
drug_trajectory = zeros(num_days*24,1);  % vector house full trajectory of EGRFi drug concentration over time
C_E_remaining = zeros(num_days,1);  % vector to house drug concentration that SHOULD remain at the end of each day (according to exponential decay equation)

for day = 1:num_days
    % Get amount of drug to add during current simulation period
%     prompt_C_E = strcat("Specify dose of EGFR inhibitor (uM) for day", num2str(day), ": ");
%     C_E_added = input(prompt_C_E);
    C_E_added = dose_by_day(day);
    initconds(3) = initconds(3) + C_E_added;

    % Numerically approximate solution to ODEs, use solutions and simulate over the period
    [time, state] = ode23(@pc9treat_sim, tspan, initconds);
   
    % Back out total population
    total_pop = state(:,1) + state(:,2);
    
    % Get EGFRi concentration that SHOULD remain at end of 24 hour period
    C_E_remaining(day) = initconds(3)*((1/2)^(24/48)); 
    
    % Add current state trajectories to the record
    population_trajectories((1+24*(day-1)):(1+24*(day)),:) = [state(:,1:2) total_pop];  % insert this simulation period's subset of the full trajectory into its appropriate hour range of total simulation
    drug_trajectory((1+24*(day-1)):(1+24*(day)),:) = state(:,3);
   
    % Update conditions for next iteration
    initconds(1) = state(end,1);
    initconds(2) = state(end,2);
    initconds(3) = state(end,3);
    
end


% 
% % Get hourly growth rate, derivative values at each hour
% state_derivs = zeros(size(state,1),2);  % initialize
% growth_rate  = zeros(size(state,1),2);
% for t = 1:size(time)
%     s0 = state(t,1:2);  % update initial condition
%     dstate_dt = pc9treat_sim(t,s0);  % get derivatives
%     state_derivs(t,:) = dstate_dt;  % store derivatives
% 
%     for col = 1:2
%         % Convert derivative (change/time) to percentage change (/time (hr))
%         growth_rate(t,col) = state_derivs(t,col)/state(t,col);
%     end
% end

% growth_rate_hr = zeros(size(state,1),2); % initialize
% for t = 1:(size(time)-1)
%     for col = 1:2
%         growth_rate_hr(t,col) = (state(t+1,col) - state(t,col))/state(t,col);
%     end
% end

% %% Plot
% time = [0:24*num_days]';  % Create vector for t axis
% 
% %%% Time Series Plot 1: Cell Populations vs Time
% figure('Name','PC9 ODE Model SIM Time Series')
% plot(time, population_trajectories, '-o')
% title('Vulnerable and Resistant Populations Over Time')
% xlabel('Time (hrs)')
% ylabel('Gross Population')
% legend('Vulnerable', 'Resistant', 'Total Population')
% % xlim([-10 106])
% ylim([-10 350])
% 
% %%% Time Series Plot 2: Drug Concentration vs Time
% figure('Name','Drug Concentration v Time')
% plot(time, drug_trajectory, '-x') 
% title('Concentration of EGFR Inhibitor Over Time')
% xlabel('Time (hrs)')
% ylabel('Drug Concentration (uM)')
% xlim([0 75])
% ylim([-1 11])
% 
% %%% Time Series Plot 3: Parameter Values vs Time
% % Obtain and plot trajectories for alpha, d_V and d_R
% alpha_trajectory = growth(drug_trajectory);
% d_V_trajectory = death_V(drug_trajectory);
% d_R_trajectory = death_R(drug_trajectory);
% trajectories = [alpha_trajectory d_V_trajectory d_R_trajectory];
% 
% % figure('Name', '/alpha, d_V, and d_R Trajectories')
% % plot(time, trajectories, '-o')
% % title('Hourly Parameter Rates During Simulation')
% % xlabel('Time (hrs)')
% % ylabel('Rate (Fraction of Relevant Population per hour)')
% % xlim([0 75])
% % legend('\alpha','d_V','d_R')
% 
% %%% Phase Plane Portrait
% % figure('Name', 'PC9 ODE Model SIM Phase Space')
% % plot(state(:,1), state(:,2))
% % title('Phase Plane Portrait')
% % xlabel('Vulnerable Population')
% % ylabel('Resistant Population')
% % xlim([-10 110])
% % ylim([-10 110])
% 
% %% Parameter Equations (copied from pc9treat_sim.m)
% % Transition (Progression G1 --> S/G2) Rate
% function alpha = growth(C_E)
%     % Hill Curve Parameters
%     a.E0 = 0.0909;%0.0225; % effect (growth rate) in absence of drug (fraction of population)
%     a.Einf = 0.015;%0.0130; % maximum possible effect (minimum growth rate) with drug (fraction of population)
%     a.EC50 = 0.015;%0.5; % concentration of drug at half max efficacy (uM)
%     a.hill = 1; % Hill slope coefficient
% 
%     alpha = a.E0 + (a.Einf - a.E0)./(1 + (a.EC50./C_E).^a.hill); % fraction of G1 population leaving per hour
% end
% 
% % Death Rates: Fraction of cells which die in a given hour
% function d_V = death_V(C_E)
%     % Hill Curve Parameters
%     v.E0 = 0.0010; % effect (death rate) in absence of drug (fraction of population)
%     v.Einf = 0.0040; % maximum possible effect (death rate) with drug (fraction of population)
%     v.EC50 = 0.25; % concentration of drug at half max efficacy (uM)
%     v.hill = 5; % Hill slope coefficient
%     
%     d_V = v.E0 + (v.Einf - v.E0)./(1 + (v.EC50./C_E).^v.hill);  % death fraction of vulnerable (G1) cells per hour
%     
% end
% function d_R = death_R(C_E)
%     % Hill Curve Parameters
%     r.E0 = 0.0010; % effect (death rate) in absence of drug (fraction of population)
%     r.Einf = 0.0040; % maximum possible effect (death rate) with drug (fraction of population)
%     r.EC50 = 0.50; % concentration of drug at half max efficacy (uM)
%     r.hill = 5; % Hill slope coefficient
%     
%     d_R = r.E0 + (r.Einf - r.E0)./(1 + (r.EC50./C_E).^r.hill);  % death fraction of resistant (S/G2) cells per hour
%     
% end
end
