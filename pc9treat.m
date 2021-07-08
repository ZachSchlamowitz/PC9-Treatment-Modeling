% Functions File for Pc9_ODE_Model
% Author: Zach Schlamowitz (7/7/21)
function state = pc9treat(t,S0)
    t_p = NaN; % hour of dosing with palbociclib; =NaN when not adding drug
    P_0 = 10;  % Initial concentration of palbociclib dose (uM)
    t_E = NaN; % hour of dosing with osimertinib; =NaN when not adding drug
    E_0 = 10;  % Initial concentration of osimertinib dose (uM)
    beta = 1/13; % Fraction of S/G2 cells which enter mitosis in a given hour

    C_p = palbo(t, t_p, P_0);
    alp = growth(C_p);
    
    C_E = osi(t, t_E, E_0);
    d_V = death_V(C_E);
    d_R = death_R(C_E);
    
    state = zeros(2,1);
    state(1) = -alp*S0(1) - d_V*S0(1) + 2*beta*S0(2);
    state(2) =  alp*S0(1) - d_R*S0(2) -   beta*S0(2);
end


% Viability Hill curves for Vulnerable, Resistant populations
function via_V = viability_V(C_E)
    % Viability Hill Curve Parameters
    v.E0 = 95;  % effect (viability) in absence of drug (percentage)
    v.Einf = 5; % maximum possible effect (viability) with drug (percentage)
    v.EC50 = 2; % effect (viability) of drug at half max efficacy concentration of drug (percentage)
    v.hill = 5; % Hill slope coefficient
    
    via_V = v.E0 + (v.Einf - v.E0)/(1 + (v.EC50/C_E)^v.hill);  % returns viability as a percentage
end
function via_R = viability_R(C_E)
    % Viability Hill Curve Parameters
    r.E0 = 95;  % effect (viability) in absence of drug (percentage)
    r.Einf = 5; % maximum possible effect (viability) with drug (percentage)
    r.EC50 = 2; % effect (viability) of drug at half max efficacy concentration of drug (percentage)
    r.hill = 5; % Hill slope coefficient
   
    % returns viability as a percentage
    via_R = r.E0 + (r.Einf - r.E0)/(1 + (2*r.EC50/C_E)^r.hill); % resistant population has twice the EC50 of vulnerable population
end

function alp = growth(C_p)
    alp = 0.85 - 0.07*C_p; % progression (G1 --> S/G2) fraction 
end
function d_V = death_V(C_E)
    d_V = (100 - viability_V(C_E))/100; % death fraction of vulnerable (G1) cells
end
function d_R = death_R(C_E)
    d_R = (100 - viability_R(C_E))/100; % death fraction of resistant (S/G2) cells
end

function C_p = palbo(t, t_p, P_0)
    if t>t_p
        C_p = P_0*(1/2)^(t/27);
    else
        C_p = 0;
    end
end
function C_E = osi(t, t_E, E_0)
    if t>t_E
        C_E = E_0*(1/2)^(t/48);
    else
        C_E = 0;
    end
end
