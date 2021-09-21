%% pc9treat_sim.m
% Functions File for Pc9_ODE_Model_1drugsim.m
% Author: Zach Schlamowitz (8/11/21)
function dstate_dt = pc9treat_sim(t,S0)%, t_E, E_0)
%     t_E = NaN; % hour of dosing with osimertinib; =NaN when not adding drug
%     t_E = input('Specify hour to dose with EGFR inhinibor (1-24): ');
%     E_0 = input('Specify dose of EGFR inhinibor (uM): ');
%     E_0 = 5;  % Initial concentration of osimertinib dose (uM)
    beta = 0.055;%1/13; % Fraction of S/G2 cells which enter mitosis in a given hour

%     C_E = osi(t, t_E, E_0);
    C_E = S0(3);
    alpha = growth(C_E);
    d_V = death_V(C_E);
    d_R = death_R(C_E);
    
    dstate_dt = zeros(2,1);
    dstate_dt(1) = -alpha*S0(1) - d_V*S0(1) + 2*beta*S0(2);
    dstate_dt(2) =  alpha*S0(1) - d_R*S0(2) -   beta*S0(2);
    dstate_dt(3) = log(1/2)*S0(3)*((1/2)^(t/48))*(1/48);
end

% Transition (Progression G1 --> S/G2) Rate
function alpha = growth(C_E)
    % Hill Curve Parameters
    a.E0 = 0.0909;%0.0225; % effect (growth rate) in absence of drug (fraction of population)
    a.Einf = 0.015;%0.0130; % maximum possible effect (minimum growth rate) with drug (fraction of population)
    a.EC50 = 0.015;%0.5; % concentration of drug at half max efficacy (uM)
    a.hill = 1; % Hill slope coefficient

    alpha = a.E0 + (a.Einf - a.E0)/(1 + (a.EC50/C_E)^a.hill); % fraction of G1 population leaving per hour
end

% Death Rates: Fraction of cells which die in a given hour
function d_V = death_V(C_E)
    % Hill Curve Parameters
    v.E0 = 0.0010; % effect (death rate) in absence of drug (fraction of population)
    v.Einf = 0.0040; % maximum possible effect (death rate) with drug (fraction of population)
    v.EC50 = 0.25; % concentration of drug at half max efficacy (uM)
    v.hill = 5; % Hill slope coefficient
    
    d_V = v.E0 + (v.Einf - v.E0)/(1 + (v.EC50/C_E)^v.hill);  % death fraction of vulnerable (G1) cells per hour
    
end
function d_R = death_R(C_E)
    % Hill Curve Parameters
    r.E0 = 0.0010; % effect (death rate) in absence of drug (fraction of population)
    r.Einf = 0.0040; % maximum possible effect (death rate) with drug (fraction of population)
    r.EC50 = 0.50; % concentration of drug at half max efficacy (uM)
    r.hill = 5; % Hill slope coefficient
    
    d_R = r.E0 + (r.Einf - r.E0)/(1 + (r.EC50/C_E)^r.hill);  % death fraction of resistant (S/G2) cells per hour
    
end

% Drug Levels
% function C_E = osi(t, t_E, E_0)
% %     if t>t_E
%         C_E = E_0*((1/2)^(t/48));
% %     else
% %         C_E = 0;
% %     end
% end
