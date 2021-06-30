% ODE Model of PC-9  Cell Response to Palbociclib and Osimertinib Treatment
% Author: Zach Schlamowitz

% Parameter values
t_p = NaN; % hour of dosing with palbociclib; =NaN when not adding drug
P_0 = 10;  % Initial concentration of palbociclib dose (uM)
t_E = NaN; % hour of dosing with osimertinib; =NaN when not adding drug
E_0 = 10;  % Initial concentration of osimertinib dose (uM)
beta = 1/13; % Fraction of S/G2 cells which enter mitosis in a given hour

% (Exaggerated) Parameter functions
syms alpha(C_p) d_V(C_E) d_R(C_E)
alpha(C_p) = 0.85 - 0.07*C_p; % progression (G1 --> S/G2) fraction 
d_V(C_E) = 100/(1.05+10*exp(-5*C_E)); % death fraction of vulnerable (G1) cells
d_R(C_E) = 50/(1+10*exp(-C_E)); % death fraction of resistant (S/G2) cells

% Main equations
syms V(t) R(t) C_p(t) C_E(t)
ode1 = diff(V,t) == -alpha*V - d_V*V + 2*beta*R;
ode2 = diff(R,t) ==  alpha*V - d_R*R -   beta*R;
% C_p(t) = piecewise(t<t_p, 0, t>t_p, P_0*(1/2)^(t-27));
C_p(t) = 0;
% C_E(t) = piecewise(t<t_E, 0, t>t_E, E_0*(1/2)^(t-48));
C_E(t) = 0;

% Specify initial conditions
cond1 = V(0) == 40; 
cond2 = R(0) == 60;

% Solve odes
odes = [ode1; ode2];
conds = [cond1; cond2];
[Vsol(t), Rsol(t)] = dsolve(odes, conds)




assert(d_V>0, 'death fraction d_V cannot be negative')
assert(d_R>0, 'death fraction d_R cannot be negative')
assert(alpha>0, 'progression fraction alpha cannot be negative')
