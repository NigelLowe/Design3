% Weights
% Design 3 - Tross
% 460368355

% weightClass inputs
% name, weight (kg), location fraction (0-1)(non dimensional), aircraft length [L](m)) 

% all equations from Nicolai unless otherwise specified

clear
clc
close all   

pl_vec = [500 3500]; % single payload value for run
en_vec = [34 24  51 40];  % single endurance value for run
% maritime field 2, airfield last 2

insideFuel = [];
params = [];
% t_vec = [];
% cg_percent_vec = [];

row = 0;
missionType = {'Maritime','Airfield'};
for mission = missionType
% basic parameters
    if strcmp(mission,'Maritime')
        albatross_parameters_maritime;
    elseif strcmp(mission,'Airfield')
        albatross_parameters_airfield;
    end
    row = row + 1;
    col = 0;
    
    for m = 1:length(pl_vec)
        pl_num = pl_vec(m);
        if strcmp(mission,'Airfield')
            m = m + 2;
        end
        en_num = en_vec(m);
        col = col + 1;

        plotOtherGraphs = 'no';
        StructureAnalysis; % runs prelim_report_code as well to get fuel used

        close all;
        clc;
        clearvars -except params row col mission missionType allWeight internalFuelWeight wingFuelWeight taper_r fused time_res pl_num pl_vec en_num en_vec c V rho_fuel insideFuel t_vec totalWeight_vec cg_percent_vec S AR b e cdo k TSFC prop_n empty_weight reach_toc cruise_alt loiter_point v_cruise v_loiter target_roc ld_climb cl_climb clmax clmin cd0 cd0c g rho0 TOW

        %internalFuelWeight = ceil(internalFuelWeight/100)*100 + 100;
        
        disp(mission);
        fprintf('%.0f kg payload | %.0f hr endurance\n', pl_num, en_num);

        % set default figure parameters
        set(groot,'defaultLineLineWidth',2.0,...
            'DefaultAxesFontSize', 20, ...
            'defaultLineMarkerSize',30,...
            'defaultAxesXGrid','on',...
            'defaultAxesYGrid','on')

        % convert variables to imperial for equations
        S_ft = S*10.7639; % ft^2


        L = 16; % m - fuselage length
        H = 2.5; % m - fueslage height
        cg_var = 0.45; % start cg position
        rearWing = cg_var;
        inFuel = 0.49; % internal fuel position

        %% Basic Weights
        % masses in lb for calculations
        MTOW = 22500*2.20462; % lb 
        M0 = 120/340; % max flight at sea level
        n = 4.5; % max load factor
        tc_ratio = 0.15;
        c_bar = 2.31; %mean(c); % mean aerodynamic chord of wing
        Lt = (1-rearWing)*L*3.28084; % ft - tail moment arm 
        sweep = 0; % deg
        tailAngle = 39; % deg
        tailRoot = 2*3.28084; % ft
        tailTip = 1.1*3.28084; % ft
        tailLength = 3.1*3.28084; % ft
        tailSpan = tailLength*cosd(tailAngle)*2; % ft
        tailArea = (tailRoot+tailTip)/2*tailLength*2; % ft^2
        tailAR = tailSpan^2/tailArea; 
        bht = tailSpan; % (ft) span of horizontal tail
        Sht = tailArea; % (ft^2) horizontal tail area
        trht = 0.15; % ft - horizontal tail thickness
        bvt = tailSpan; % vertical tail height
        Svt = tailArea; % ft^2 - vertical tail area
        Sr = 1.55*0.67 * 3.28084^2; % ft^2 - rudder area (Sr/Svt=0.3 in unknown)
        ARvt = bvt^2/Svt; % vertical tail aspect ratio
        taper_rvt = tailTip/tailRoot; % vertical tail taper ratio
        sweepvt = 10; % deg - vertical tail sweep

        Kinl = 1.25; % for inlets in fuselage
        q = 0.5*rho0*V^2*0.02088547; % lb/ft^2 - max dynamic pressure

        Mland = MTOW/1.5; % lb - max landing weight
        Ln = 1.2*39.3701; % in - nose landing gear
        Lm = 1.45*39.3701; % in - main landing gear weight 
        %fuselageWeight1 = 10.43*Kinl^1.42*(q*10^-2)^0.245*(MTOW*10^-3)^0.98*(L/H)^0.71; % USAF/Commercial - Nicolai
        fuselageWeight = 11.03*Kinl^1.23*(q*10^-2)^0.245*(MTOW*10^-3)^0.98*(L/H)^0.61; % USN - Nicolai
        wingWeight = 1.75 * 0.00428*S_ft^0.48*(AR*M0^0.43*(MTOW*n)^0.84*taper_r^0.14)/((100*tc_ratio)^0.76*(cosd(sweep))^1.54); % Subsonic Aircraft - Nicolai % multiply by 1.75 since the C-2 Grumman has a similar wingspan and folding wings with outer weight 3000lb.
        % https://www.tested.com/tech/568755-airplane-origami-how-folding-wings-work/     ^
        vertTailWeight = 0.0034*((MTOW*n)^0.813*Sht^0.584*(bht/trht)^0.033*(c_bar*Lt)^0.28)^0.915; % Nicolai - scale since have V tail
        horiTailWeight = 0.19*((MTOW*n)^0.363*Svt^1.089*M0^0.601*Lt^(-0.726)*(1+Sr/Svt)^0.217*ARvt^0.337*(1+taper_rvt)^0.363*(cosd(sweepvt)^(-0.484)))^1.014; % Nicolai - scale since have v tail
        tailWeight = vertTailWeight*sind(tailAngle)^2+horiTailWeight*cosd(tailAngle)^2; % lb - scale since using V tail https://vtechworks.lib.vt.edu/bitstream/handle/10919/26482/Dissertation3a.pdf?sequence=1&isAllowed=y
        frontGearWeight = 0.125*(Mland*n)^0.566*(Ln/12)^0.845; % Nicolai according to Gudmundsson
        rearGearWeight = 2*0.054*(Mland*n)^0.684*(Lm/12)^0.601; % Nicolai according to Gudmundsson

        basic(1) = weightClass(              'Fuselage',  fuselageWeight,      0.5, L);
        basic(2) = weightClass(             'Main Wing',      wingWeight, rearWing, L); % multiplied by 1.75 for folding mechanism
        basic(3) = weightClass(                  'Tail',      tailWeight,     0.96, L);
        basic(4) = weightClass(    'Front Landing Gear', frontGearWeight,      0.2, L);
        basic(5) = weightClass(     'Rear Landing Gear',  rearGearWeight,      0.7, L);
        basic.lb2kg; % convert masses from lb to kg
        [basicWeight, basicMoment] = basic.totalWM;
        fprintf('total basic:      %.0f kg | %.0f kgm\n', basicWeight, basicMoment);

        %% Propulsion Weights
        % masses in lb for calculations
        Ni = 2; % number of inlets
        Ld = 1; % ft - subsonic duct length
        Ai = 0.5; % ft^2 - area per inlet
        P2 = 10.49707; % lb/in^2 - max static pressure at engine compressor face
        Kgeo = 1; % duct shape factor (1 for round duct)
        Km = 1; % duct material factor (1 for M < 1)
        %Lr = 0.6; % ft - ramp length forward of throat per inlet
        engineWeight = 4189; % lb - TP-400
        ne = 1; % number of engines
        np = 1; % number of propellers
        nb = 4; % number of blades per propeller
        dp = 3*3.28084; % ft - propeller diameter
        hp = 11000; % shaft horsepower
        %Kte = 1; % temperature correction factor (1 for Max nach number < 1)
        duct = 0.32*Ni*Ld*Ai^0.65*P2*0.6; % duct support structure (internal only)
        internalDuct = 1.735*(Ni*Ld*Ai*0.5*P2*Kgeo*Km)^0.7331; % duct structure from inlet lip to engine compressor face (internal engine installations only)
        totalDuct = duct + internalDuct;
        wingFuel = sum(wingFuelWeight);
        totalFuelWeight = wingFuel+internalFuelWeight;
        maxFuelWeight = 18341; % max fuel - so the below vaues are the same for all cases
        totalFuelGallons = maxFuelWeight/rho_fuel*1000*0.214172; % Imperial gallon (should be US (0.26) but dont want a bigger value)
        bladderCells = 23.1*(totalFuelGallons*10^-2)^0.758; % non-self sealing bladder cells
        cellSupports = 7.91*(totalFuelGallons*10^-2)^0.854; % fuel system bladder cell backing and supports
        totalBladder = bladderCells + cellSupports;
        dumpDrain = 7.38*(totalFuelGallons*10^-2)^0.458; % dump and drain system
        cgFuelControl = 7.91*(totalFuelGallons*10^-2)^0.854; % cg control system (transfer pumps and monitor)
        pneumaticTP = 12.05*(ne*engineWeight*10^-3)^1.458;
        propellerControls = 0.322*nb^0.589*(np*dp*hp*10^-3)^1.178; % propeller controls - turboprop engine

        % rotor estimation
        % https://rotorcraft.arc.nasa.gov/Publications/files/NDARCTheory_v1_6_938.pdf
        Nrotor = 2;
        Nblade = 2;
        shaftWallT = 1.2*MTOW/1000;
        liftOffset = 0.15;
        rotorRadius = 6.75*3.28084;
        separationFraction = 1/(6.75*2);
        tipClearance = separationFraction*0.1;
        Vtip = 984.252;
        bladeTR = 0.9;
        t2r = 0.12*(0.8+0.2*bladeTR)/(0.5+0.5*bladeTR); 
        
        w_blade = 4*113.3*2.20462;
        %w_blade = Nrotor*0.00008377*shaftWallT*liftOffset*rotorRadius^3/(2*(separationFraction - tipClearance)*t2r^2);
        w_hub = Nrotor*(0.17153*shaftWallT*rotorRadius*Nblade + 0.000010534*(w_blade/Nrotor)*Vtip^2*t2r/rotorRadius);
        w_shaft = Nrotor*0.081304*shaftWallT*liftOffset*rotorRadius^2*2*separationFraction/t2r;
        
        prop(1)  = weightClass(               'Engine',    engineWeight,  0.19, L);
        prop(2)  = weightClass(     'Propeller Blades',     100*2.20462, 0.015, L); % 374.1795 lb - using equation
        prop(3)  = weightClass('Propeller Drive Shaft',       0*2.20462,     0, L);

        prop(4)  = weightClass(              'Duct',         totalDuct,   0.05, L);
        prop(5)  = weightClass(     'Pneumatic (TP)',       pneumaticTP,    0.2, L);
        prop(6)  = weightClass( 'Propeller Controls', propellerControls,    0.2, L);
        prop(7)  = weightClass('Fuel System Bladder',      totalBladder, inFuel, L);
        prop(8)  = weightClass(   'Fuel Dump System',         dumpDrain,   0.67, L);
        prop(9) = weightClass(      'Fuel CG Control',    cgFuelControl,   0.75, L);

        if strcmp(mission, 'Maritime')  %cg_var
            rotorPos = 0.45;
            prop(10)  = weightClass( 'Rotor Blades',     w_blade,      rotorPos, L); 
            prop(11)  = weightClass(    'Rotor Hub',       w_hub,      rotorPos, L);
            prop(12)  = weightClass(  'Rotor Shaft',     w_shaft,      rotorPos, L);
            prop(13)  = weightClass('Rotor gearbox', 300*2.20462, rotorPos*0.85, L);
        end
        prop.lb2kg; % convert masses from lb to kg
        [propWeight, propMoment] = prop.totalWM;
        fprintf('total propulsion: %.0f kg | %.0f kgm\n', propWeight, propMoment);

        %% Payload Weights
        % masses in kg
        % fixed equipment weights (check if greater than 500 kg internal payload
        payload(1) = weightClass(         'MX-20',                          90,      0.9, L);
        payload(2) = weightClass(      'Lynx SAR',                          38,     0.85, L);
        payload(3) = weightClass('Other Internal PL', 500-payload(1:2).totalWM,      0.8, L); % the rest of the 500kg internal payload
        payload(4) = weightClass(      'External PL',               pl_num-500, rearWing, L);
        [payloadWeight, payloadMoment] = payload.totalWM;
        fprintf('total payload:     %.0f kg | %.0f kgm\n', payloadWeight,payloadMoment);

        %% Fuel Weight
        % masses in kg %inFuel
        inFuel = 0.49;
        fuelStart(1) = weightClass(      'Wing Fuel',     wingFuelWeight, rearWing, L);
        fuelStart(2) = weightClass(  'Internal Fuel', internalFuelWeight,   inFuel, L);
        fuel(1) = weightClass(    'Wing Fuel Flight',     wingFuelWeight, rearWing, L);
        fuel(2) = weightClass('Internal Fuel Flight', internalFuelWeight,   inFuel, L);
        [fuelWeight, fuelMoment] = fuelStart.totalWM;
        fprintf('total fuel:      %.0f kg | %.0f kgm\n\n',fuelWeight,fuelMoment);    

        %%
        % empty (no fuel)
        totalWeight0 = basicWeight + propWeight;
        totalMoment0 = basicMoment + propMoment;
        cg_empty0 = (totalMoment0/totalWeight0)/L*100;
        
        totalWeight = basicWeight + propWeight + payloadWeight;
        totalMoment = basicMoment + propMoment + payloadMoment;
        cg_empty = totalMoment/totalWeight;

        % with fuel for this mission
        totalWeight(1) = totalWeight + fuelWeight;
        totalMoment = totalMoment + fuelMoment;
        fprintf('total empty: %.0f kg\n',basicWeight + propWeight);
        fprintf('total weight: %.0f kg\n',totalWeight(1));
        fprintf('total moment: %.0f kgm\n\n',totalMoment);

        cg(1) = totalMoment/totalWeight(1);
        cg_percent = cg(1)/L*100;
        x_ac = 7.1;
        fprintf('cg empty: %.3f (%.3f %%)(%.3f %%)\n',cg_empty,cg_empty/L*100,(cg_empty-x_ac)/c_bar);
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
        
        insideFuel = [insideFuel; fuelStart(1).weight fuelStart(2).weight fuel(1).weight internalFuelWeight wingFuel totalFuelWeight totalFuelWeight/rho_fuel]; % inside fuel (kg),total wing fuel (kg), total fuel (kg), total fuel (m^3)
        
        params = [params; cg_empty0, cg_empty/L*100, cg(1)/L*100, totalWeight(1)]; % oew_cg (%), zfw_cg (%), start_cg, TOW, 
        
        t_vec{row,col} = (0:length(fused))*time_res;
        cg_percent_vec{row,col} = cg/L*100; %(cg - x_ac)/c_bar; %
        allWeight{row,col} = totalWeight;
        
    end
end

%%
arraySize = size(cg_percent_vec,1);
xLimits = [43.5 47];
yLimits1 = [0 55];
yLimits2 = [5000 25000];
figure(1)
subplot(1,2,1)
% yyaxis left
for j = 1:arraySize
%     yyaxis left
    plot(cg_percent_vec{1,j},t_vec{1,j},'linewidth',3)
    legend_labels{j} = sprintf('%.0f kg.', pl_vec(j));
    
%     yyaxis right
%     cg_percent_length = (cg_percent_vec{1,j}*c_bar + x_ac)/L; 
%     set(gca,'YTick',cg_percent_length, 'YLim', [min(cg_percent_length), max(cg_percent_length)]) % #important
%     ylabel('cg location (% aft of nose)');
    
    hold on
end
% yyaxis left
legend(legend_labels,'location','NE')
ylabel('endurance (hr)')
xlabel('cg location (% aircraft length)')
xlim(xLimits)
ylim(yLimits1)
title(missionType{1})


subplot(1,2,2)
for j = 1:arraySize %arraySize/2+1:arraySize
    plot(cg_percent_vec{2,j},t_vec{2,j},'linewidth',3)
    hold on
end
legend(legend_labels,'location','NE')
ylabel('endurance (hr)')
xlabel('cg location (% aircraft length)')
xlim(xLimits)
ylim(yLimits1)
title(missionType{2})



figure(2)
subplot(1,2,1)
for j = 1:arraySize
    plot(cg_percent_vec{1,j},allWeight{1,j},'linewidth',3)
    legend_labels{j} = sprintf('%.0f kg.', pl_vec(j)); 
    hold on
end
legend(legend_labels,'location','SE')
ylabel('weight (kg)')
xlabel('cg location (% aircraft length)')
xlim(xLimits)
ylim(yLimits2)
title(missionType{1})


subplot(1,2,2)
for j = 1:arraySize
    plot(cg_percent_vec{2,j},allWeight{2,j},'linewidth',3)
    hold on
end
legend(legend_labels,'location','SE')
ylabel('weight (kg)')
xlabel('cg location (% aircraft length)')
xlim(xLimits)
ylim(yLimits2)
title(missionType{2})



% Print order of components
allArray = [basic prop payload fuelStart];
totalSize = allArray.length;
allWeights(1:totalSize) = allArray;

keys = cell(1,totalSize);
[values,valuesW] = deal(zeros(1,totalSize));
for i = 1:totalSize
    keys{i} = allWeights(i).name;
    values(i) = allWeights(i).location/L;
    valuesW(i) = allWeights(i).weight;
end
[InOrderLocation, sortIdx] = sort(values);
InOrderNames = keys(sortIdx);
InOrderWeights = valuesW(sortIdx);
front2BackLocations = [InOrderNames; num2cell(InOrderLocation); num2cell(InOrderWeights)];
disp('       Object name     | location (n.d.) | weight (kg)')
disp(front2BackLocations');

