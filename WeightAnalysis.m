% Weights
% Design 3 - Tross
% 460368355

% weightClass inputs
% name, weight (kg), location fraction (0-1)(non dimensional), aircraft length [L](m)) 

% all equations from Nicolai unless otherwise specified

clear

pl_vec = [3500 500]; % single payload value for run
en_vec = [30   45];  % single endurance value for run

t_vec = [];
totalWeight_vec = [];
cg_percent_vec = [];

for k = 1:length(pl_vec)
    pl_num = pl_vec(k);
    en_num = en_vec(k);
    
    plotOtherGraphs = 'no';
    StructureAnalysis; % runs prelim_report_code as well to get fuel used

    close all;
    clc;
    %clearvars -except internalFuelWeight wingFuelWeight beamWeight taper_r fused time_res pl_num en_num c V rho_fuel

    fprintf('%.0f kg payload | %.0f hr endurance\n', pl_num, en_num);
    
    % set default figure parameters
    set(groot,'defaultLineLineWidth',2.0,...
        'DefaultAxesFontSize', 20, ...
        'defaultLineMarkerSize',30,...
        'defaultAxesXGrid','on',...
        'defaultAxesYGrid','on')

    % basic parameters
    albatross_parameters_maritime;
    %albatross_parameters_airfield;

    % convert variables to imperial for equations
    S = S*10.7639; % ft^2


    L = 16; % m - fuselage length
    H = 2.5; % m - fueslage height
    rearWing = 0.45;
    cg_var = 0.45; % start cg position
    inFuel_loc = 0.45; % internal fuel position

    %% Basic Weights
    % masses in lb for calculations
    MTOW = 20000*2.20462; % lb 
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
    tailAngle = 30; % deg
    Mland = MTOW/2; % lb - max landing weight
    Ln = 39.3701; % in - nose landing gear
    Lm = 39.3701; % in - main landing gear weight 
    %fuselageWeight = 10.43*Kinl^1.42*(q*10^-2)^0.245*(MTOW*10^-3)^0.98*(L/H)^0.71; % USAF/Commercial - Nicolai
    fuselageWeight = 11.03*Kinl^1.23*(q*10^-2)^0.245*(MTOW*10^-3)^0.98*(L/H)^0.61; % USN - Nicolai
    wingWeight = 0.00428*S^0.48*(AR*M0^0.43*(MTOW*n)^0.84*taper_r^0.14)/((100*tc_ratio)^0.76*(cosd(sweep))^1.54); % Subsonic Aircraft - Nicolai
    vertTailWeight = 0.0034*((MTOW*n)^0.813*Sht^0.584*(bht/trht)^0.033*(c_bar*Lt)^0.28)^0.915; % Nicolai - scale since have V tail
    horiTailWeight = 0.19*((MTOW*n)^0.363*Svt^1.089*M0^0.601*Lt^(-0.726)*(1+Sr/Svt)^0.217*ARvt^0.337*(1+taper_rvt)^0.363*(cosd(sweepvt)^(-0.484)))^1.014; % Nicolai - scale since have v tail
    tailWeight = vertTailWeight*sind(tailAngle)^2+horiTailWeight*cosd(tailAngle)^2; % lb - scale since using V tail https://vtechworks.lib.vt.edu/bitstream/handle/10919/26482/Dissertation3a.pdf?sequence=1&isAllowed=y
    frontGearWeight = 0.125*(Mland*n)^0.566*(Ln/12)^0.845; % Nicolai according to Gudmundsson
    rearGearWeight = 0.054*(Mland*n)^0.684*(Lm/12)^0.601; % Nicolai according to Gudmundsson

    basic(1) = weightClass(          'Fuselage',  fuselageWeight,      0.5, L);
    basic(2) = weightClass(         'Main Wing',      wingWeight, rearWing, L);
    basic(3) = weightClass(              'Tail',      tailWeight,     0.95, L);
    basic(4) = weightClass('Front Landing Gear', frontGearWeight,      0.2, L);
    basic(5) = weightClass( 'Rear Landing Gear',  rearGearWeight,      0.8, L);
    basic.lb2kg; % convert masses from lb to kg
    [basicWeight, basicMoment] = basic.totalWM;
    fprintf('total basic:      %.0f kg | %.0f kgm\n', basicWeight, basicMoment);

    %% Propulsion Weights
    % masses in lb for calculations
    Ni = 2; % number of inlets
    Ld = 1; % ft - subsonic duct length
    Ai = 0.5; % ft^2 - area per inlet
    P2 = 1; % lb/in^2 - max static pressure at engine compressor face
    Kgeo = 1; % duct shape factor (1 for round duct)
    Km = 1; % duct material factor (1 for M < 1)
    Lr = 0.6; % ft - ramp length forward of throat per inlet
    engineWeight = 4189; % lb
    ne = 1; % number of engines
    np = 1; % number of propellers
    nb = 3; % number of blades per propeller
    dp = 3*3.28084; % ft - propeller diameter
    hp = 11000; % shaft horsepower
    Kte = 1; % temperature correction factor (1 for Max nach number < 1)
    duct = 0.32*Ni*Ld*Ai^0.65*P2*0.6; % duct support structure (internal only)
    internalDuct = 1.735*(Ni*Ld*Ai*0.5*P2*Kgeo*Km)^0.7331; % duct structure from inlet lip to engine compressor face (internal engine installations only)
    totalDuct = duct + internalDuct;
    engineControl = 4.079*(Ni*Lr*Ai^0.5*Kte)^1.201;
    totalFuelGallons = (sum(wingFuelWeight)+internalFuelWeight)/rho_fuel*1000*0.214172; % Imperial gallon (should be US (0.26) but dont want a bigger value)
    bladderCells = 23.1*(totalFuelGallons*10^-2)^0.758; % non-self sealing bladder cells
    cellSupports = 7.91*(totalFuelGallons*10^-2)^0.854; % fuel system bladder cell backing and supports
    totalBladder = bladderCells + cellSupports;
    dumpDrain = 7.38*(totalFuelGallons*10^-2)^0.458; % dump and drain system
    cgFuelControl = 7.91*(totalFuelGallons*10^-2)^0.854; % cg control system (transfer pumps and monitor)
    pneumaticTP = 12.05*(ne*engineWeight*10^-3)^1.458;
    propellerControls = 0.322*nb^0.589*(np*dp*hp*10^-3)^1.178; % propeller controls - turboprop engine

    prop(1)  = weightClass(               'Engine',    engineWeight,    0.3, L);
    prop(2)  = weightClass(     'Propeller Blades',      20*2.20462,      0, L);
    prop(3)  = weightClass(    'Propeller Gearbox',     100*2.20462,   0.05, L);
    prop(4)  = weightClass('Propeller Drive Shaft',      30*2.20462,   0.18, L);
    prop(5)  = weightClass(         'Rotor Blades', 4*113.3*2.20462, cg_var, L); % need to adjust time calculation for storing rotor 
    prop(6)  = weightClass(        'Rotor Gearbox',     750*2.20462, cg_var*1.25, L);
    prop(7)  = weightClass(          'Rotor Shaft',      120*2.20462,cg_var, L);

    prop(8)  = weightClass(                 'Duct',         totalDuct,          0, L);
    prop(9)  = weightClass(       'Engine Control',     engineControl,        0.7, L);
    prop(10) = weightClass(       'Pneumatic (TP)',       pneumaticTP,       0.15, L);
    prop(11) = weightClass(   'Propeller Controls', propellerControls,       0.15, L);
    prop(12) = weightClass(  'Fuel System Bladder',      totalBladder, inFuel_loc, L);
    prop(13) = weightClass(     'Fuel Dump System',         dumpDrain,   rearWing, L);
    prop(14) = weightClass(      'Fuel CG Control',     cgFuelControl,       0.7, L);

    prop.lb2kg; % convert masses from lb to kg
    [propWeight, propMoment] = prop.totalWM;
    fprintf('total propulsion: %.0f kg | %.0f kgm\n', propWeight, propMoment);

    %% Payload Weights
    % masses in kg
    % fixed equipment weights (check if greater than 500 kg internal payload
    payload(1) = weightClass(         'MX-20',                          90,      0.8, L);
    payload(2) = weightClass(      'Lynx SAR',                          57,      0.7, L);
    payload(3) = weightClass('Other Internal PL', 500-payload(1:2).totalWM,      0.9, L); % the rest of the 500kg internal payload
    payload(4) = weightClass(      'External PL',               pl_num-500, rearWing, L);
    [payloadWeight, payloadMoment] = payload.totalWM;
    fprintf('total payload:     %.0f kg | %.0f kgm\n', payloadWeight,payloadMoment);

    %% Fuel Weight
    % masses in kg
    fuelStart(1) = weightClass(      'Wing Fuel',     wingFuelWeight,   rearWing, L);
    fuelStart(2) = weightClass(  'Internal Fuel', internalFuelWeight, inFuel_loc, L);
    fuel(1) = weightClass(    'Wing Fuel Flight',     wingFuelWeight,   rearWing, L);
    fuel(2) = weightClass('Internal Fuel Flight', internalFuelWeight, inFuel_loc, L);
    [fuelWeight, fuelMoment] = fuelStart.totalWM;
    fprintf('total fuel:      %.0f kg | %.0f kgm\n\n',fuelWeight,fuelMoment);    
%     fuelStart(1) = weightClass(      'Wing Fuel',     0,   rearWing, L);
%     fuelStart(2) = weightClass(  'Internal Fuel',    0, inFuel_loc, L);
%     fuel(1) = weightClass(    'Wing Fuel Flight',     0,   rearWing, L);
%     fuel(2) = weightClass('Internal Fuel Flight', 0, inFuel_loc, L);
%     [fuelWeight, fuelMoment] = fuelStart.totalWM;

    %%
    totalWeight(1) = basicWeight + propWeight + payloadWeight + fuelWeight;
    totalMoment = basicMoment + propMoment + payloadMoment + fuelMoment;
    fprintf('total empty: %.0f kg\n',basicWeight + propWeight);
    fprintf('total weight: %.0f kg\n',totalWeight(1));
    fprintf('total moment: %.0f kgm\n',totalMoment);

    cg(1) = totalMoment/totalWeight(1);
    cg_percent = cg(1)/L*100;
    fprintf('cg start: %.3f m (%.3f %%)\n\n',cg(1),cg_percent);

    for i = 1:length(fused)
        fuelU = fused(i);

        if fuel(2).weight > 0 % use internal fuel first 
            fuel(2).updateWeight(fuel(2).weight - fuelU);
        else % use wing  fuel after internal finished
            fuel(1).updateWeight(fuel(1).weight - fuelU);
        end

        if i == round(length(fused)/2)
            payload(4).updateWeight(0);
            [payloadWeight, payloadMoment] = payload.totalWM;
        end

        [fuelWeight, fuelMoment] = totalWM(fuel);  
        totalWeight(i+1) = basicWeight + propWeight + payloadWeight + fuelWeight;
        totalMoment = basicMoment + propMoment + payloadMoment + fuelMoment;
        cg(i+1) = totalMoment/totalWeight(i+1);

    end

    t_vec = [t_vec; (0:length(fused))*time_res];
    totalWeight_vec = [totalWeight_vec; totalWeight];
    cg_percent_vec = [cg_percent_vec; cg/L*100];
end
    


% figure(1)
% for j = 1:size(cg_percent_vec,1)
%     plot(cg_percent_vec(j,:),totalWeight_vec(j,:),'linewidth',3)
%     legend_labels{j} = sprintf('%.0f hrs.', en_vec(j));
%     hold on
% end
% legend(legend_labels,'location','NW')
% % plot(cg_percent_vec,totalWeight_vec);
% xlabel('cg location aft of nose (%)')
% ylabel('Aircraft Weight (kg)');

figure(1)
% t = (0:length(fused))*time_res;
% plot(t,cg_percent_vec);
for j = 1:size(cg_percent_vec,1)
    plot(t_vec(j,:),cg_percent_vec(j,:),'linewidth',3)
    legend_labels{j} = sprintf('%.0f kg.', pl_vec(j));
    hold on
end
legend(legend_labels,'location','NW')
xlabel('endurance (hr)')
ylabel('cg location aft of nose (%)')


% Print order of components
allArray = [basic prop payload fuelStart];
totalSize = allArray.length;
allWeights(1:totalSize) = allArray;

keys = cell(1,totalSize);
values = zeros(1,totalSize);
for i = 1:totalSize
    keys{i} = allWeights(i).name;
    values(i) = allWeights(i).location/L;
    valuesW(i) = allWeights(i).weight;
end
[InOrderLocation, sortIdx] = sort(values);
InOrderNames = keys(sortIdx);
InOrderWeights = valuesW(sortIdx);
front2BackLocations = [InOrderNames; num2cell(InOrderLocation); num2cell(InOrderWeights)];
disp('Object location (non-dimensional)')
disp(front2BackLocations');