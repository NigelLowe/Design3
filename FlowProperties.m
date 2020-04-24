% function to calculate air density and dynamic pressure from altitude and
% airspeed
function [T,P,rho,Mach,q_bar] = FlowProperties(h,V)
% constants
gamma = 1.4;
Trate=0.0065; 
R=287.05; 
% base values at 0m
T0=288.15;  
P0=101325;


% gravity
g = 9.81;

T=T0-Trate*h;
P=P0*(1-Trate*h/T0)^(g/(R*Trate));
rho=P/(R*T);

% speed of sound
a = sqrt(gamma*R*T);
% Mach number
Mach = V/a;
% dynamic pressure    
q_bar = 0.5*rho*V^2;
end

