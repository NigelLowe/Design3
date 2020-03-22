clear
clc
clf
close all;

h_cruise = 40000*0.3048;
% Wing characteristics
S = 120; % Wing area (m^2)
WS = 35; % Wingspan (m)
AR = (WS^2)/S; % Aspect ratio ((span^2)/(wing area))
Clmax = 1.8;
Clmin = 0.1;
Cd0 = 0.06; %0.026; % Drag coefficient at 0 angle of attack
e = 0.75; % Oswald efficiency
K = 1/(pi()*AR*e);

% Constants
g = 9.81; % acceleration due to gravity (m/s^2)
mu = 0.02; % Friction of runway
R = 286.9; % Gas constant for air (J/kg.K)
rho0 = 1.225; % Density at sea level (kg/m^3)
Laps = 0.0065; % (K/m)Lapse rate

% Weight characteristics
Wmtow = 25000; % MTOW mass (lb)
Wempty = 8700; % empty mass (lb)
W = Wmtow*g; % MTOW (N)

% fuel = 128; % (gallons) fuel capacity
% rho_w = 1000; % density of water
% fsg = 0.718; % fuel specific gravity
% fuel = fuel*(3.78541e-3)*rho_w*fsg; % (kg) fuel capacity
% fuel_lb = fuel/2.20462;
% fuel_flow = 19.2; % gph
% fuel_flow = fsg*fuel_flow*3.785/3600; % kg/s
% 
% Fvariable = 0.15*fuel; % variable fuel reserve
% Freserve = fuel_flow*45*60; % required reserve fuel
% Wfinal = W - (fuel - Fvariable - Freserve)*g; % final allowable weight
Wfinal = 8700;

% Propeller characteristics
Power0 = 11000*745.7; % shaft power at sea level
% TSFC = 0.39;%/3600; (N/hr / N)Thrust Specific Fuel Consumption
TSFC = 0.85*18*1e-6 * g*3600;   %kg/(s.N)

% Sea Level characteristics
T0 = 15 + 273.15;
P0 = 101.3e3;

% % Read data from blade analysis program
% fileID = fopen('propeller.txt','r');
% formatSpec = '%d %f';
% sizeA = [8 Inf];
% A = fscanf(fileID,formatSpec,sizeA);
% A = A';
% fclose('all');

% Altitude and Velocity range
h_max = 50000; % ft
h = 0:250:h_max; % altitudes (ft)
v_max = 300*0.514444; % m/s
V = 0:v_max; % (m/s) 
V_knots = V*1.94384; % (knots)

% J = V/(nD)
% Work out velocity from J value
% V2 = A(:,2)*n*d;
% efficiency = A(:,8); 
% eff = max(efficiency) - 0.1;
eff = 0.8;

theta = 0;
count = 0;


for i=1:length(h)
    alt = h(i)*0.3048; % (m) altitudes 
    temp = T0 - alt*Laps; % (K) temperature at these altitudes
    press=P0*(temp/T0)^(g/Laps/R); % (Pa) pressures at these altitudes
    rho = press/(R*temp); % (kg/m^3) density at these altitudes 
    P = Power0*rho/rho0; % (W) power adjusted for altitude
 
    % Stall boundary
    vstall(i)=sqrt(W/Clmax/0.5/rho/S);
    
    for j = 1:length(V)
        Vel = V(j);
        
        T = P*eff/Vel;       % (N) thrust
        Cl = W/(0.5*rho*Vel^2*S); % lift coefficient
        Cd = Cd0 + K*(Cl-Clmin)^2;        % drag coefficient
        D = 0.5*Cd*rho*Vel^2*S;   % (N) drag

        % Specific Excess Power
        Ps(i,j) = Vel*(T - D)/W; % W would change throughout the flight meanign so would L and D. Need to add a time loop to calculate
        
        % Range
        Range(i,j) = (Vel/(TSFC/3600))*(Cl/Cd)*(log(W/Wfinal));% (m) /1000000;        
        
        time(i,j) = Range(i,j)/Vel;
        
        % Specific energy
        he(i,j) = alt + 0.5*(V(j)^2/g);
        
        % calculate average climb angle
        if h(i) <= h_cruise && V(j) >= vstall(i)
            angle_climb(i,j) = (T - D)/W;
        end
    end        
end

%%

% Contour Specific Excess Power
figure(1);
subplot(2,2,1)
vals = 0:50;
[c1,h1]=contour(V_knots,h,Ps,vals,'k');
clabel(c1,h1, 'labelspacing',1e10);
title('Specific Excess Power Contours');
xlabel('Velocity (knots)');
ylabel('Altitude (ft)');
hold on;
plot (vstall*1.94384,h,'-r');
legend('Ps','stall boundary');

% Contour Range
% figure(2)
subplot(2,2,2)
vals2 = 0:1000:500000;
[c2,h2]=contour(V_knots,h,Range/1000,vals2,'k'); % Multiply by 0.621371 - convert to miles
clabel(c2,h2, 'labelspacing',1e10);
xlabel('Velocity (knots)')
ylabel('Altitude (ft)')
title('Range Contours');
hold on
[c3,h3] = contour(V_knots,h,Ps,[1 1000]);
clabel(c3,h3, 'labelspacing',1e10);
plot(vstall*1.94384,h,'-r')
legend('Range contours (km)','Ps=1','stall boundary')

% Contour Time
% figure(3)
subplot(2,2,3)
vals3 = 0:10:100;
[c4,h4]=contour(V_knots,h,time/3600,vals3,'k'); % Multiply by 0.621371 - convert to miles
clabel(c4,h4, 'labelspacing',1e10);
xlabel('Velocity (knots)')
ylabel('Altitude (ft)')
title('Time Contours');
hold on
[c5,h5] = contour(V_knots,h,Ps,[1 1000]);
clabel(c5,h5, 'labelspacing',1e10);
plot(vstall*1.94384,h,'-r')
legend('Time contours (hr)','Ps=1','stall boundary')