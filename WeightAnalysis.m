% Weights
% Design 3 - Tross
% 460368355

clear

pl_num = 500; %3500;
en_num = 40; %20; 

plotOtherGraphs = 'no';
StructureAnalysis; % runs prelim_report_code as well to get fuel used

close all;
clc;
clearvars -except internalFuelWeight wingFuelWeight beamWeight taper_r fused time_res pl_num en_num c V

% set default figure parameters
set(groot,'defaultLineLineWidth',2.0,...
    'DefaultAxesFontSize', 20, ...
    'defaultLineMarkerSize',30,...
    'defaultAxesXGrid','on',...
    'defaultAxesYGrid','on')

% basic parameters
albatross_parameters;

% class inputs
% name, weight (kg), location fraction (0-1)(non dimensional), aircraft length [L](m)) 


% convert variables to imperial for equations
S = S*10.7639; % ft^2


L = 16; % m - fuselage length
H = 2.5; % m - fueslage height
rearWing = 0.45;
cg_var = 0.54;
inFuel_loc = 0.5;

%% Basic Weights
% masses in lb for calculations
MTOW = 15000*2.20462; % lb 
M0 = 120/340; % max flight at sea level
n = 4.5; % max load factor
tc_ratio = 0.15;
sweep = 0; % deg
Sv = 7; % V-Tail area
ARt = 3; % Tail aspect ratio
Sht = 1*7*sind(30)*2; % (ft^2) horizontal tail area
bht = 7*sind(30)*2; % (ft) span of horizontal tail
trht = 0.5; % ft - horizontal tail thickness
c_bar = mean(c); % mean aerodynamic chord of wing
Lt = (1-rearWing)*L*3.28084; % ft - tail moment arm 
Svt = 1*7*cosd(30)*2; % ft^2 - vertical tail area
bvt = 7*cosd(30); % vertical tail height
Sr = 0.3*Svt; % ft^2 - rudder area (Sr/Svt=0.3 in unknown)
ARvt = bvt^2/Svt; % vertical tail aspect ratio
taper_rvt = 0.8; % vertical tail taper ratio
sweepvt = 10; % deg - vertical tail sweep
Kinl = 1.25; % for inlets in fuselage
q = 0.5*rho0*V^2*0.02088547; % lb/ft^2 - max dynamic pressure

%fuselageWeight = 10.43*Kinl^1.42*(q*10^-2)^0.245*(MTOW*10^-3)^0.98*(L/H)^0.71; % USAF/Commercial - Nicolai
fuselageWeight = 11.03*Kinl^1.23*(q*10^-2)^0.245*(MTOW*10^-3)^0.98*(L/H)^0.61; % USN - Nicolai
wingWeight = 0.00428*S^0.48*(AR*M0^0.43*(MTOW*n)^0.84*taper_r^0.14)/((100*tc_ratio)^0.76*(cosd(sweep))^1.54); % Subsonic Aircraft - Nicolai
vertTailWeight = 0.0034*((MTOW*n)^0.813*Sht^0.584*(bht/trht)^0.033*(c_bar*Lt)^0.28)^0.915; % Nicolai - scale since have V tail
horiTailWeight = 0.19*((MTOW*n)^0.363*Svt^1.089*M0^0.601*Lt^(-0.726)*(1+Sr/Svt)^0.217*ARvt^0.337*(1+taper_rvt)^0.363*(cosd(sweepvt)^(-0.484)))^1.014; % Nicolai - scale since have v tail
tailWeight = 0.75*(vertTailWeight+horiTailWeight); % lb - scale since using V tail
%landingGearWeight = 62.21*(MTOW*10^-3)^0.84; % USAF/Commercial - Nicolai
landingGearWeight = 129.1*(MTOW*10^-3)^0.66; % USN - Nicolai
frontGearWeight = 0.3*landingGearWeight;
rearGearWeight = 0.7*landingGearWeight;

basic(1) = weightClass(          'Fuselage',  fuselageWeight,      0.5, L);
basic(2) = weightClass(         'Main Wing',      wingWeight, rearWing, L);
basic(3) = weightClass(              'Tail',      tailWeight,     0.95, L);
basic(4) = weightClass('Front Landing Gear', frontGearWeight,      0.2, L);
basic(5) = weightClass( 'Rear Landing Gear',  rearGearWeight,      0.7, L);
basic.lb2kg; % convert masses from lb to kg
[basicWeight, basicMoment] = basic.totalWM;
fprintf('total basic:      %.0f kg | %.0f kgm\n', basicWeight, basicMoment);

%% Propulsion Weights
% still need to add extra stuff for propulsion

prop(1) = weightClass(               'Engine',    1900,    0.3, L);
prop(2) = weightClass(     'Propeller Blades',      20,      0, L);
prop(3) = weightClass(    'Propeller Gearbox',     200,   0.05, L);
prop(4) = weightClass('Propeller Drive Shaft',     100,   0.17, L);
prop(5) = weightClass(         'Rotor Blades', 4*113.3, cg_var, L); % need to adjust time calculation for storing rotor 
prop(6) = weightClass(        'Rotor Gearbox',    1000, cg_var, L);
prop(7) = weightClass(          'Rotor Shaft',     200, cg_var, L);
[propWeight, propMoment] = prop.totalWM;
fprintf('total propulsion: %.0f kg | %.0f kgm\n', propWeight, propMoment);

%% Payload Weights
payload(1) = weightClass(         'MX-20',                   90,      0.8, L);
payload(2) = weightClass(      'Lynx SAR',                   57,      0.1, L);
payload(3) = weightClass('Other Internal', 500-totalWM(payload),     0.25, L);
payload(4) = weightClass(      'External',           pl_num-500, rearWing, L);
[payloadWeight, payloadMoment] = totalWM(payload);
fprintf('total payload:    %.0f kg | %.0f kgm\n', payloadWeight,payloadMoment);

%% Fuel Weight
fuelStart(1) = weightClass(      'Wing Fuel',     wingFuelWeight,   rearWing, L);
fuelStart(2) = weightClass(  'Internal Fuel', internalFuelWeight, inFuel_loc, L);
fuel(1) = weightClass(    'Wing Fuel Flight',     wingFuelWeight,   rearWing, L); %fuelStart(1);
fuel(2) = weightClass('Internal Fuel Flight', internalFuelWeight, inFuel_loc, L); %fuelStart(2);
[fuelWeight, fuelMoment] = totalWM(fuelStart);
fprintf('total fuel:       %.0f kg | %.0f kgm\n\n',fuelWeight,fuelMoment);    

%%
totalWeight(1) = basicWeight + propWeight + payloadWeight + fuelWeight;
totalMoment = basicMoment + propMoment + payloadMoment + fuelMoment;
fprintf('total empty:    %.0f kg\n',basicWeight + propWeight);
fprintf('total weight: %.0f kg\n',totalWeight(1));
fprintf('total moment: %.0f kgm\n',totalMoment);

cg(1) = totalMoment/totalWeight(1);
cg_percent = cg(1)/L*100;
fprintf('cg: %.3f m\n',cg(1));
fprintf('cg: %.3f %%\n',cg_percent);
    
for i = 1:length(fused)
    fuelU = fused(i);
   
    if fuel(2).weight > 0 % use internal fuel first 
        fuel(2).updateWeight(fuel(2).weight - fuelU);
    else % use wing  fuel equally
        fuel(1).updateWeight(fuel(1).weight - fuelU);
    end
    
    if i == round(length(fused)/2)
        payload(4).weight = 0;
    end
        
    [fuelWeight, fuelMoment] = totalWM(fuel);  
    totalWeight(i+1) = basicWeight + propWeight + payloadWeight + fuelWeight;
    totalMoment = basicMoment + propMoment + payloadMoment + fuelMoment;
    cg(i+1) = totalMoment/totalWeight(i+1);
    
end

cg_percent = cg/L*100;


figure(1)
plot(cg_percent,totalWeight);
xlabel('cg location aft of nose (%)')
ylabel('Aircraft Weight (kg)');

figure(2)
t = 0:time_res:en_num+time_res;
plot(t,cg_percent);
xlabel('endurance (hr)')
ylabel('cg location aft of nose (%)')
