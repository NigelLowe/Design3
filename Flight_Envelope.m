% Flight envelope
% 460389730

% clear command window and workspace
clc;
clear;
clf;

% Input aerodynamic parameters.
MTOW = 1435;            % lbs
S = 191.5;              % ft^2 wing area
b = 57.4;               % wing span in feet    
rho = 0.0023769;        % slugs/ft^3 density factor at service ceiling*density atsea level
V_stalldown = 38;       %  knots
V_stallup = 35;         %  knots
W_L = MTOW/S;           % lb/ft^2 wing loading
C = S/b;                % mean chord, feet
AR = b^2/S;             % wing aspect ratio
g = 32.2;               % gravity in imperial ft/s^2
U_1 = 50;               % 50 ft/s gust cruise
U_2 = 25;               % 25 ft/s gust dive

%% Calculations

% Max lift at stall
CLmax_Down = MTOW/(0.5*rho*S*(V_stalldown^2));
CLmax_Up = MTOW/(0.5*rho*S*(V_stallup^2));

% Designed cruise speed.
V_C = 92;

% Limiting dive speed
V_D = 135;

% Lift curve slope
a = (2*pi)/(1+2/AR);
ug = (2*W_L)/(rho*C*g*a);
K = (0.88 * ug)/(5.3 +ug);

% max loading factor
nmax_positive = 4.2; %2.1 + (24000/(MTOW + 10000));
nmax_negative = -1*0.4 * nmax_positive;

% maximum load factor due to cruise gust or dive gust
gust_cruise_pos = 1+((K*U_1*V_C*a)/(498*W_L));
gust_cruise_neg = 1-((K*U_1*V_C*a)/(498*W_L));
gust_dive_pos = 1+((K*U_2*V_D*a)/(498*W_L));
gust_dive_neg = 1-((K*U_2*V_D*a)/(498*W_L));

% Max and Min Manouverability speeds.
Va_Up = sqrt((nmax_positive*MTOW)/(0.5*rho*CLmax_Up*S)); 
Va_Down = sqrt((-1*nmax_negative*MTOW)/(0.5*rho*CLmax_Down*S)); 

%% Array for the cruise manourevility
n_p= 0:0.001:nmax_positive;
Va_p = sqrt((n_p(:)*MTOW)./(0.5*rho*CLmax_Up*S));

% plot arrays
figure(1);
hold on;
plot(Va_p,n_p,'b');
title('Flight Envelope');
xlabel('Airspeed (knots)');
ylabel('Load Factor n');

text(10,4,'Gust Limits','Color','red','FontSize',10)
text(10,3.5,'Manoeuvre Limits','Color','blue','FontSize',10)


%% Array for the dive manourevility
n_n = 0:-0.001:nmax_negative;
Va_n = sqrt(abs((n_n(:)*MTOW)./(0.5*rho*CLmax_Down*S)));
plot(Va_n,n_n,'b');

%% Top line from A to D
vAD = [Va_Up V_D];
nAD = [nmax_positive nmax_positive];
plot(vAD, nAD,'b');

%% Back line from D to E
vDE = [V_D V_D];
nDE = [nmax_positive 0];
plot(vDE, nDE,'b');

%% Bottom line from G to F
vGD = [Va_Down V_C];
nGD = [nmax_negative nmax_negative];
plot(vGD, nGD,'b');

%% Oblique back line F to E
vFE = [V_C V_D];
nFE = [nmax_negative 0];
plot(vFE, nFE,'b')

%% V Stall lines.
Vsx = [V_stallup V_stallup];
Vsy = [1 0];
plot(Vsx, Vsy, 'b');

Vsx1 = [V_stalldown V_stalldown];
Vsy1 = [0 -1];
plot(Vsx1, Vsy1, 'b');


%% Plot gust limits

% gust for positve at 50knots.
vGustCruise1 = [0 V_C];
nGustCruise1 = [1 gust_cruise_pos];
plot(vGustCruise1, nGustCruise1,'r');

% gust for negative at 50 knots.
vGustCruise2 = [0 V_C];
nGustCruise2 = [1 gust_cruise_neg];
plot(vGustCruise2, nGustCruise2,'r');

% gust for positve at 25knots.
vGustDive1 = [0 V_D];
nGustDive1 = [1 gust_dive_pos];
plot(vGustDive1, nGustDive1,'r');

% gust for negative at 25knots.
vGustDive2 = [0 V_D];
nGustDive2 = [1 gust_dive_neg];
plot(vGustDive2, nGustDive2,'r');

%% Gust limit filling lines
% top oblique
v_gCD = [V_C V_D];
n_gCD = [gust_cruise_pos gust_dive_pos];
plot(v_gCD, n_gCD,'r');

% bottom oblique
v_gFD = [V_C V_D];
n_gFD = [gust_cruise_neg gust_dive_neg];
plot(v_gFD, n_gFD,'r');

% Back line
v_gFD = [V_D V_D];
n_gFD = [gust_dive_neg 0];
plot(v_gFD, n_gFD,'r');

%% Axis lines
x1 = [V_C V_C];
y1 = [gust_cruise_pos gust_cruise_neg];
plot(x1,y1,'k');

x2 = [0 V_D];
y2 = [0 0];
plot(x2,y2,'k');

x = [0 V_D];
y = [1 1];
plot(x,y,'k');

axis([0 140 -3.0482 5.0482]);
box on
