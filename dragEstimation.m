% Drag profile estimation
% Author: Jack Knight
% Date: 06/04/2020

% Clear the workspace
clear;
clc;

% Get the aircraft parameters
albatross_parameters_maritime;

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

% Choose the wing area as the reference area
Sref = S;

%% Get the wetted area for all area perpendicular to flow
t = 0.152*c;    % Max thickness of airfoil
L_D = 15;       % Best lift to drag ratio

Mach = v_mean/a;    % Mach number during cruise

% Force terms
f_tc = 1 + 2.7*(t/c) + 100*(t/c)^4;
f_m = 1 - 0.08*Mach^1.45;

% Fuselage
L_fuse = 16;
D_fuse = 2;
f_ld = 1 + 60/((L_fuse/D_fuse)^2) + 0.0025*(L_fuse/D_fuse);
S_fuse = pi * (D_fuse/2)^2;
Cdo_fuse = Cf * f_ld * f_m * S_fuse/Sref;
fprintf('\nFuselage: Swet = %.4f m^2, Cdo = %.5f\n', S_fuse, Cdo_fuse)

% Wing 
S_wing = b*t;
Cdmin_airfoil = 0.01;
Cdo_wing = Cf * f_tc * f_m * S_wing/Sref * (Cdmin_airfoil/0.004)^0.4;
fprintf('Wing:     Swet = %.4f m^2, Cdo = %.5f\n', S_wing, Cdo_wing)

% Tailplane
b_tail = 4.6;
t_tail = 0.152 * 0.5*(2.5+1.1);
S_tail = b_tail * t_tail;
Cdo_tail = Cf * f_tc * f_m * S_tail/Sref * (Cdmin_airfoil/0.004)^0.4;
fprintf('Tail:     Swet = %.4f m^2, Cdo = %.5f\n', S_tail, Cdo_tail)

% Rotor


% Payload - External Pods
no_pods = 4;
D_pod = 0.18;
L_pod = 1.8;
S_pod = pi * (D_pod/2)^2;
f_ld_pod = 1 + 60/((L_pod/D_pod)^2) + 0.0025*(L_pod/D_pod);
Cdo_pod = Cf * f_ld_pod * f_m * S_pod/Sref;
fprintf('Pod:      Swet = %.4f m^2, Cdo = %.5f\n', S_fuse, no_pods*Cdo_pod)

%Interference Factor
Q_nacelle = 1.5*3;
Q_wing = 1;
Q_fuselage = 1.04;
Q_fuselage = 1;

Q_total = Q_nacelle+Q_wing+Q_fuselage+Q_fuselage;

% Total skin friction Drag
Cdf = Q_total * (Cdo_fuse + Cdo_wing + Cdo_tail + no_pods*Cdo_pod);
fprintf('\nSkin Friction Drag coef, Cf = %.4f\n', Cdf)

%Cdf = Cf * Swet / Sref;
Cdo = 1.2 * Cdf;
fprintf('Min Drag Coef, Cdo = %.4f\n', Cdo)