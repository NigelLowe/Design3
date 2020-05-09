% Flight envelope
% 460389730

% certify under STANAG Subpart C 321-

% clear command window and workspace
clc;
clear;
clf;

% set default figure parameters
        set(groot,'defaultLineLineWidth',2.0,...
            'DefaultAxesFontSize', 20, ...
            'defaultLineMarkerSize',30,...
            'defaultAxesXGrid','on',...
            'defaultAxesYGrid','on')
        
%albatross_parameters_maritime;
albatross_parameters_airfield;

% Input aerodynamic parameters.
g = 9.81;                 % gravity in SI
MTOW = 6300*g;           % N --------------------------------------------- from design min to design max weight 
rho0 = 1.225;             % kg/m3 - SL density 
W_L = MTOW/S;             % kg/m^2 wing loading

alt = 0000; % ft
if alt < 20000
    rho = rho0;
    U_1 = 15.2;               % m/s - 50 ft/s gust cruise
    U_2 = 7.6;                % m/s - 25 ft/s gust dive
else
    U_1 = 15.2 - 19/75000 * (alt - 20000);
    U_2 = 7.6 - 19/150000 * (alt - 20000);
    rho = rho0*3.468/14.696; % 35000ft factor    -------------------------- needs to be done at range of operating altitude
end

%% Calculations

% Max lift at stall
CLmax_Up = 2; %MTOW/(0.5*rho*S*(V_stallup^2));
CLmax_Down = 1.5; %MTOW/(0.5*rho*S*(V_stalldown^2));
V_stallup = sqrt(MTOW/(0.5*rho*CLmax_Up*S)); %52;   % m/s
V_stalldown = sqrt(MTOW/(0.5*rho*CLmax_Down*S)); %55; % m/s

% Designed cruise speed.
V_C = 200*0.514444; % knots

% Limiting dive speed
V_D = 320*0.514444; % knots

% max loading factor
nmax_positive = min(3.8, 2.1 + (10900/(22500*g/g + 4536)));
nmax_negative = -0.4 * nmax_positive;

% Lift curve slope
a = (2*pi)/(1+2/AR);
ug = (2*W_L)/(rho*c*g*a);
K = (0.88 * ug)/(5.3 + ug);

% maximum load factor due to cruise gust or dive gust
gust_cruise_pos = 1+((K*rho0*U_1*V_C*a)/(2*W_L));
gust_cruise_neg = 1-((K*rho0*U_1*V_C*a)/(2*W_L));
gust_dive_pos   = 1+((K*rho0*U_2*V_D*a)/(2*W_L));
gust_dive_neg   = 1-((K*rho0*U_2*V_D*a)/(2*W_L));

% Max and Min Manouverability speeds.
Va_Up = sqrt((nmax_positive*MTOW)/(0.5*rho*CLmax_Up*S)); 
Va_Down = sqrt((-1*nmax_negative*MTOW)/(0.5*rho*CLmax_Down*S)); 

%% Array for the cruise manourevility
n_p = 0:0.001:nmax_positive;
Va_p = sqrt((n_p(:)*MTOW)./(0.5*rho*CLmax_Up*S));

% Convert velocities to knots for plotting
Va_p = Va_p*1.94384;
Va_Up = Va_Up*1.94384;
Va_Down = Va_Down*1.94384;
V_C = V_C*1.94384;
V_D = V_D*1.94384;
V_stallup = V_stallup*1.94384;
V_stalldown = V_stalldown*1.94384;

% plot arrays
figure(1);
hold on;
plot(Va_p,n_p,'b');
%title('Flight Envelope');
xlabel('Airspeed (knots)');
ylabel('Load Factor n');
text(10,3.5,'Gust Limits','Color','red','FontSize',15)
text(10,3,'Manoeuvre Limits','Color','blue','FontSize',15)


%% Array for the dive manourevility
n_n = 0:-0.001:nmax_negative;
Va_n = sqrt(abs((n_n(:)*MTOW)./(0.5*rho*CLmax_Down*S)));
Va_n = Va_n*1.94384; % convert to knots
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

%axis([0 V_D*1.1 -3.0482 5.0482]);
%box on
ylim([-2 4])