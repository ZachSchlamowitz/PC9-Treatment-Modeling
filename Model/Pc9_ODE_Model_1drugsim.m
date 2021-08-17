%% Pc9_ODE_Model_1drugsim.m
% 1 Drug ODE Model of PC-9 Cell Response to treatment with varying doses of EGFR inhibition
% Author: Zach Schlamowitz (8/11/21)

%% Main equations

% tspan = [0 240];
tspan = [0:72];
initconds = [40; 60];
[time, state] = ode23(@pc9treat_sim, tspan, initconds);
state(:,3) = state(:,1) + state(:,2);  % total population

% Get hourly growth rate, derivative values at each hour
state_derivs = zeros(size(state,1),2);  % initialize
growth_rate  = zeros(size(state,1),2);
for t = 1:size(time)
    s0 = state(t,1:2);  % update initial condition
    dstate_dt = pc9treat_sim(t,s0);  % get derivatives
    state_derivs(t,:) = dstate_dt;  % store derivatives

    for col = 1:2
        % Convert derivative (change/time) to percentage change (/time (hr))
        growth_rate(t,col) = state_derivs(t,col)/state(t,col);
    end
end

% growth_rate_hr = zeros(size(state,1),2); % initialize
% for t = 1:(size(time)-1)
%     for col = 1:2
%         growth_rate_hr(t,col) = (state(t+1,col) - state(t,col))/state(t,col);
%     end
% end

% Plot
figure('Name','PC9 ODE Model SIM Time Series')
plot(time, state, '-o')
title('Vulnerable and Resistant Populations Over Time')
xlabel('Time (hrs)')
ylabel('Gross Population')
legend('Vulnerable', 'Resistant', 'Total Population')
% xlim([-10 106])
ylim([-10 350])

% Plot Phase Plane
% figure('Name','PC9 ODE Model SIM Phase Space')
% plot(state(:,1), state(:,2))
% title('Phase Plane Portrait')
% xlabel('Vulnerable Population')
% ylabel('Resistant Population')
% xlim([-10 110])
% ylim([-10 110])

