% ODE Model of PC-9  Cell Response to Palbociclib and Osimertinib Treatment
% Author: Zach Schlamowitz

%% Parameter values
% t_p = NaN; % hour of dosing with palbociclib; =NaN when not adding drug
% P_0 = 10;  % Initial concentration of palbociclib dose (uM)
% t_E = NaN; % hour of dosing with osimertinib; =NaN when not adding drug
% E_0 = 10;  % Initial concentration of osimertinib dose (uM)
% beta = 1/13; % Fraction of S/G2 cells which enter mitosis in a given hour
% 
% % %% (Exaggerated) Parameter functions
% % % Structs for viabilty Hill curve parameters
% v.E0 = 95;  % effect (viability) in absence of drug (percentage)
% v.Einf = 5; % maximum possible effect (viability) with drug (percentage)
% v.EC50 = 2; % effect (viability) of drug at half max efficacy concentration of drug (percentage)
% v.hill = 5; % Hill slope coefficient
% r = v;  % NOTE: here MATLAB creates separate copy of struct rather than an alias
% 
% % Viability Hill curves for Vulnerable, Resistant populations
% syms viability_V(C_E) viability_R(C_E)
% viability_V(C_E) = v.E0 + (v.Einf - v.E0)/(1 + (v.EC50/C_E)^v.hill);
% viability_R(C_E) = v.E0 + (v.Einf - v.E0)/(1 + (2*v.EC50/C_E)^v.hill); % resistant population has twice the EC50 of vulnerable population
% 
% syms alpha(C_p) d_V(C_E) d_R(C_E)
% alpha(C_p) = 0.85 - 0.07*C_p; % progression (G1 --> S/G2) fraction 
% d_V(C_E) = 100 - viability_V(C_E); % death fraction of vulnerable (G1) cells
% d_R(C_E) = 100 - viability_R(C_E); % death fraction of resistant (S/G2) cells
% 
% 
% %% FIRST ATTEMPT (FAILED/INCOMPLETE):
% syms C_p(t) C_E(t)
% % ode1 = diff(V,t) == -alpha*V - d_V*V + 2*beta*R;
% % ode2 = diff(R,t) ==  alpha*V - d_R*R -   beta*R;
% % % C_p(t) = piecewise(t<t_p, 0, t>t_p, P_0*(1/2)^(t-27));
% C_p(t) = 0;
% % % C_E(t) = piecewise(t<t_E, 0, t>t_E, E_0*(1/2)^(t-48));
% C_E(t) = 0;
% % 
% % % Specify initial conditions
% % cond1 = V(0) == 40; 
% % cond2 = R(0) == 60;
% % 
% % % Solve odes
% % odes = [ode1; ode2];
% % conds = [cond1; cond2];
% % [Vsol(t), Rsol(t)] = dsolve(odes, conds)  %z: Andrew suggests using ode23 rather than dsolve. This uses numerical approach.
% % 
% 
% 
% 
% % assert(d_V>0, 'death fraction d_V cannot be negative')
% % assert(d_R>0, 'death fraction d_R cannot be negative')
% % assert(alpha>0, 'progression fraction alpha cannot be negative')

%% Main equations

tspan = [0 72];
tspan = [0:72];
initconds = [40; 60];
[time, state] = ode23(@pc9treat, tspan, initconds);
state(:,3) = state(:,1) + state(:,2);  % total population

% Get hourly growth rate, derivative values at each hour
state_derivs = zeros(size(state,1),2);  % initialize
growth_rate  = zeros(size(state,1),2);
for t = 1:size(time)
    s0 = state(t,1:2);  % update initial condition
    dstate_dt = pc9treat(t,s0);  % get derivatives
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
figure('Name','PC9 ODE Model Time Series')
plot(time, state, '-o')
title('Vulnerable and Resistant Populations Over Time')
xlabel('Time (hrs)')
ylabel('Gross Population')
legend('Vulnerable', 'Resistant', 'Total Population')
xlim([-10 106])
ylim([-10 350])

% Plot Phase Plane
figure('Name','PC9 ODE Model Phase Space')
plot(state(:,1), state(:,2))
title('Phase Plane Portrait')
xlabel('Vulnerable Population')
ylabel('Resistant Population')
% xlim([-10 110])
% ylim([-10 110])


% function state = pc9treat(t,S0)
%     state = zeros(2,1);
%     state(1) = -alpha*S0(1) - d_V*S0(1) + 2*beta*S0(2);
%     state(2) =  alpha*S0(1) - d_R*S0(2) -   beta*S0(2);
% end
