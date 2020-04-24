% Drag profile estimation
% Author: Jack Knight
% Date: 06/04/2020

% Clear the workspace
if ~exist('plotOtherGraphs','var') % if statement for this file to be used in other functions
    clear
    
    % Call parameters
    albatross_parameters_maritime
    %albatross_parameters_airfield
end 
clc;

%% Calculate atmospheric Properties

% Standard atmosphere - assume linearly decreasing temperature
L = 0.0065;     % Lapse rate (K/m)
g = 9.81;       % Gravitational acceleration (m/s^2)
gamma = 1.4;    % Ratio of specific heats
R = 287;        % Gas Constant for air

% Conditions at Sea Level
T_sl = 15 + 273.15; % Temperature (K)
P_sl = 101300;      % Pressure (Pa)

% Air Density at sea level (kg/m^3)
rho_sl = P_sl / (R*T_sl); 

% Determine conditions at Altitude
T = T_sl - (cruise_alt * L);            
P = P_sl * (T/T_sl)^(g/(L*R));  
rho = P / (R*T); 
a = sqrt(gamma*R*T);

% Viscosity
mu = 1.7894e-5;
nu = mu/rho;

%% Simple approximation - Reynolds number

% Flight speed - mean cruise speed
v_mean = 220/1.944;

Re = c*v_mean/nu;

Cf_lam = 1.327 / (sqrt(Re));
Cf_turb = 0.455 / (log10(Re)^2.58);

% Check if flow is laminar or turbulent
if Re > 4000
    fprintf('Turbulent flow, Re = %.2f\n', Re)
    Cf = Cf_turb;
elseif Re < 2100
    fprintf('Laminar flow, Re = %.2f\n', Re)
    Cf = Cf_lam;
else
    fprintf('Transitional flow, Re = %.2f\n', Re)
    Cf = (Cf_lam + Cf_turb)/2;
end

fprintf('Skin Friction Coefficient, Cf = %.4f\n', Cf)



%% Get the wetted area for all area perpendicular to flow

% Which config?
config = 1;
t = 0.152*c;    % Max thickness of airfoil
L_D = 15;       % Best lift to drag ratio

% Choose the wing area as the reference area
Sref = c*b/1.45;

Mach = v_mean/a;    % Mach number during cruise

% Config chech
if config == 0
    no_pods = 0;
    no_rotor = 0;
else
    no_pods = 4;
    no_rotor = 2;
end

% Force terms
f_tc = 1 + 2.7*(t/c) + 100*(t/c)^4;
f_m = 1 - 0.08*Mach^1.45;

% Fuselage
L_fuse = 16;
D_fuse = 2;
f_ld = 1 + 60/((L_fuse/D_fuse)^2) + 0.0025*(L_fuse/D_fuse);
%S_fuse = pi * (D_fuse/2)^2;
S_fuse = pi*D_fuse*L_fuse*(0.5+0.135*(4/L_fuse))^2/3*(1.015+0.3/(L_fuse/D_fuse)^1.5);
Cdo_fuse = Cf * f_ld * f_m * S_fuse/Sref;
fprintf('\nFuselage: Swet = %6.3f m^2, Cdo = %.5f\n', S_fuse, Cdo_fuse)

% Wing 
lambda = 1.2/2.7; %Taper Ratio
%S_wing = b*c;
S_wing = 2*b*c*(1+0.25*(t/c)*((1+t/c*lambda)/(1+lambda)));
Cdmin_airfoil = 0.01;
Cdo_wing = Cf * f_tc * f_m * S_wing/Sref * (Cdmin_airfoil/0.004)^0.4;
fprintf('Wing:     Swet = %6.3f m^2, Cdo = %.5f\n', S_wing, Cdo_wing)

% Tailplane
b_tail = 4.6;
c_tail = 0.5*(2.5+1.1);
t_tail = 0.152 * c_tail;
lambda_tail = 1.1/2.5;
%S_tail = b_tail * c_tail;
S_tail = 2*b_tail*c_tail*(1+0.25*(t/c)*((1+t/c*lambda_tail)/(1+lambda_tail)));
Cdo_tail = Cf * f_tc * f_m * S_tail/Sref * (Cdmin_airfoil/0.004)^0.4;
fprintf('Tail:     Swet = %6.3f m^2, Cdo = %.5f\n', S_tail, Cdo_tail)

% Rotor
b_rotor = 6.75*2;
c_rotor = 0.7;
t_rotor = 0.12 * c_rotor;
%S_rotor = b_rotor * c_rotor;
S_rotor = 2*b_rotor*c_rotor*(1+0.25*(t_rotor/c_rotor)*((1+t_rotor/c_rotor*1)/(1+1)));
Q_rotor = 1.2;
Cdo_rotor = no_rotor * Q_rotor * (Cf * f_tc * f_m * S_rotor/Sref * (0.02/0.004)^0.4);
fprintf('Rotor:    Swet = %6.3f m^2, Cdo = %.5f\n', S_rotor, Cdo_rotor)

% Payload - External Pods
D_pod = 0.4;
L_pod = 1.8;
%S_pod = pi * (D_pod/2)^2;
S_pod = pi*D_pod*L_pod*(0.5+0.135*(4/L_pod))^2/3*(1.015+0.3/(L_pod/D_pod)^1.5);
f_ld_pod = 1 + 60/((L_pod/D_pod)^2) + 0.0025*(L_pod/D_pod);
Q_pod = 1.5;
Cdo_pod = Cf * Q_pod * f_ld_pod * f_m * S_pod/Sref;
fprintf('Pod:      Swet = %6.3f m^2, Cdo = %.5f\n', S_pod, Cdo_pod)

fprintf('\nAirfield Config:')
% Total skin friction Drag
Cdf = 1 * (Cdo_fuse + Cdo_wing + Cdo_tail + no_pods*Cdo_pod);
fprintf('\nSkin Friction Drag coef, Cf = %.4f\n', Cdf)

%Cdf = Cf * Swet / Sref;
Cdo = 1.2 * Cdf;    % Extra 20% Accounts for pressure drag
Cdo = 1.05 * Cdo;   % Extra 5% Accounts for Interference effects
fprintf('Min Drag Coef, Cdo = %.4f\n', Cdo)

%% Maritime config
fprintf('\nMaritime Config:')
% Total skin friction Drag
Cdf = 1 * (Cdo_fuse + Cdo_wing + Cdo_tail + Cdo_rotor + no_pods*Cdo_pod);
fprintf('\nSkin Friction Drag coef, Cf = %.4f\n', Cdf)

%Cdf = Cf * Swet / Sref;
Cdo = 1.2 * Cdf;    % Extra 20% Accounts for pressure drag
Cdo = 1.05 * Cdo;   % Extra 5% Accounts for Interference effects
fprintf('Min Drag Coef, Cdo = %.4f\n', Cdo)

%% Angus equation

%Interference Factor
Q_nacelle = 1.5*3;
Q_wing = 1;
Q_fuselage = 1.04;
Q_fuselage = 1;

Q_total = Q_nacelle+Q_wing+Q_fuselage+Q_fuselage;

% Form Factor
xi = 0.4; %position of max thickness
sweep = 1; % No sweep, cos(x) = 1
FF = (1+(0.6/xi)*(t/c)+100*(t/c)^4)*(1.34*Mach^0.18*(sweep)^0.28);

Swet_total = S_fuse + S_wing + S_tail + no_pods*S_pod;
CD0_new = Cf * FF * Q_total * (Swet_total/Sref);
%fprintf('Angus Drag: Cdo = %.4f\n', CD0_new)
