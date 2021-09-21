%% Pc9_ODE_Model_1drugsim.m
% 1 Drug ODE Model of PC-9 Cell Response to treatment with varying doses of EGFR inhibition
% Author: Zach Schlamowitz (8/11/21)

%% Main equations

% User Initializations:
num_days = 3;
% C_E = 0;
initconds = [40; 60; 0];  % number of cells in [G1 S/G2] and EGFRi concentration [EGFRi]

% Internal Initializations:
% tspan = [0 72];
tspan = [0:24];
full_state = zeros(num_days*24,3);
full_drug_state = zeros(num_days*24,1);
C_E_remaining = zeros(num_days,1);

for day = 1:num_days
    % Get hour of dosing and add dose amount to current cumulative dose
%     prompt_t_E = strcat("Specify hour to dose with EGFR inhibitor (1-24) for Day ", num2str(day), ": "); % FLAG: need to figure out what happens when t<t_E / is it a problem
%     t_E = input(prompt_t_E); % +1? To account for fact that matlab vars are not zero based?
%     C_E = C_E + input("Specify dose of EGFR inhibitor (uM): ");
    prompt_C_E = strcat("Specify dose of EGFR inhibitor (uM) for day", num2str(day), ": ");
    C_E_added = input(prompt_C_E);
    initconds(3) = initconds(3) + C_E_added;
%     [time, state] = ode23(@(tspan, initconds) pc9treat_sim(tspan, initconds, t_E, C_E), tspan, initconds);
    [time, state] = ode23(@pc9treat_sim, tspan, initconds);
    % Back out total population
%     state(:,3) = state(:,1) + state(:,2);
    total_pop = state(:,1) + state(:,2);
    % Get remaining EGFRi concentration at end of 24 hour period
    C_E_remaining(day) = initconds(3)*((1/2)^(24/48));  % FLAG: need to figure out how to get C_E from one day to carry over to next without overwriting or setting C_E to 0
    
    % Add current state trajectories to the record
    full_state((1+24*(day-1)):(1+24*(day)),:) = [state(:,1:2) total_pop];
    full_drug_state((1+24*(day-1)):(1+24*(day)),:) = state(:,3);
    % Update conditions for next iteration
    initconds(1) = state(end,1);
    initconds(2) = state(end,2);
    initconds(3) = state(end,3);
    
end

time = [0:24*num_days]';

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

% Plot
figure('Name','PC9 ODE Model SIM Time Series')
plot(time, full_state, '-o')
title('Vulnerable and Resistant Populations Over Time')
xlabel('Time (hrs)')
ylabel('Gross Population')
legend('Vulnerable', 'Resistant', 'Total Population')
% xlim([-10 106])
ylim([-10 350])

figure('Name','Drug Concentration v Time')
plot(time, full_drug_state, '-x') 
title('Concentration of EGFR Inhibitor Over Time')
xlabel('Time (hrs)')
ylabel('Drug Concentration (uM)')
xlim([0 75])
ylim([-1 11])


% Plot Phase Plane
% figure('Name', 'PC9 ODE Model SIM Phase Space')
% plot(state(:,1), state(:,2))
% title('Phase Plane Portrait')
% xlabel('Vulnerable Population')
% ylabel('Resistant Population')
% xlim([-10 110])
% ylim([-10 110])

