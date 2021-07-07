% ODE Practice Script
% Zach Schlamowitz, 07/06/21
%% Ex.1: Using ode23/ode45 to solve simple equation y'=2t (from documentation 'ode23')
% tspan = [0 5];
% y0 = 0;
% func = @(t,y) 2*t;
% [t,y] = ode23(func, tspan, y0);
% 
% plot(t,y, '-o')

%% Ex.2: Van der Pol, x'' - mu(1-x^2)x' + x = 0, 
%                i.e. x' = y, y' = mu(1-x^2)y - x 
% (from documentation 'ode23')

% % Define ode solver parameters
% tspan = [0 20];  % time points
% y0 = [2; 0];  % initial condition
% [t,y] = ode45(@vanderpol, tspan, y0);
% 
% % Plot time series
% plot(t, y(:,1), '-o', t, y(:,2), '-o')
% title('Solution of van der Pol Equation (\mu = 1) with ode45')
% xlabel('Time t')
% ylabel('Solution y')
% legend('y_1', 'y_2')
% 
% % Function defs
% function dydt = vanderpol(t,y)
%     dydt = [y(2); 1*(1-y(1)^2)*y(2) - y(1)];
% end

%% Ex.3: Lotka-Volterra (from documentation 'Solve Predator-Prey Equations')
% x' = x - axy
% y' = -y + Bxy
% 
% t0 = 0;
% tf = 15;
% initvals = [20; 20];
% [time, sol] = ode23(@lotka, [t0 tf], initvals);
% 
% % Time Series Plot
% plot(time, sol)
% title('Predator, Prey Populations Over Time')
% xlabel('t')
% ylabel('Population')
% legend('Prey', 'Predators', 'Location', 'North')
% 
% % Phase Plane Plot
% plot(sol(:,1), sol(:,2))
% title('Phase Space')
% xlabel('Prey Population')
% ylabel('Predator Population')
% 
% % Add in version with ode45 as well
% [T,S] = ode45(@lotka, [t0 tf], initvals);
% plot(sol(:,1), sol(:,2), '-', S(:,1), S(:,2), '-');
% title('Comparing Phase Plane Portraits')
% legend('ode23', 'ode45')
% 
% function eqs = lotka(t,initvals)
%     a = 0.01;
%     B = 0.02;
%     eqs = zeros(2,1);
%     eqs(1) = initvals(1) - a*initvals(1)*initvals(2);
%     eqs(2) = -initvals(2) + B*initvals(1)*initvals(2);
% end

%% Ex.4: My Own Creation, an autonomous, nonlinear system
% x' = 3x + x^2 + xy
% y' = xy + 2y^2

% tspan = [0 25];
% initconds = [3; 0.5];
% [time, sol] = ode23(@equ, tspan, initconds);
% 
% % Plot Time Series
% plot(time, sol, '-o')
% legend('People', 'Martians')
% title('Example ODE System')
% xlabel('time')
% ylabel('Population')
% 
% % Plot Phase Plane
% plot(sol(:,1), sol(:,2))
% title('Phase Plane')
% xlabel('People')
% ylabel('Martians')
% 
% function z = equ(t,z0)
%     z = zeros(2,1);
%     z(1) = 3*z0(1) - z0(1)^2 - z0(1)*z0(2);
%     z(2) = -z0(1)*z0(2) - 2*z0(2)^2;
% end
