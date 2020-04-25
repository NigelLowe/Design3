% Structural Analysis
% 460368355

if ~exist('plotOtherGraphs','var') % if statement for this file use in other functions
clear
    
    % general parameters
    albatross_parameters_maritime;
    albatross_parameters_airfield;
end
clc
close all

% set default figure parameters
set(groot,'defaultLineLineWidth',2.0,...
    'DefaultAxesFontSize', 20, ...
    'defaultLineMarkerSize',30,...
    'defaultAxesXGrid','on',...
    'defaultAxesYGrid','on')

%% general parameters
if ~exist('plotOtherGraphs','var') % if statement for this file use in other functions
    pl_num = 500;
    en_num = 51;
end
plotOtherGraphs = 'no';
prelim_report_code;
dragEstimation;
clc
close all

% for taper ratio input
%taper_r = 0.45; % from graph in lecture 5 slide 19. Minimum drag at 0 sweep
% predator - taper_r = 0.25
% reaper - taper_r = 0.8
% global hawk - taper_r = 0.3
%{ 
https://www.semanticscholar.org/paper/A-Genetic-Algorithm-Incorporating-Design-Choice-for-Mull/0742bbd9eb7ea5651c01241ae2357bde480a4625?tab=abstract&citingPapersSort=is-influential&citingPapersLimit=10&citingPapersOffset=0&year%5B0%5D=&year%5B1%5D=&citedPapersSort=is-influential&citedPapersLimit=10&citedPapersOffset=10
%}

cr = 2*b/(AR*(1+taper_r)); % m - chord at tip
ct = cr * taper_r; % m - chord at root


% equations of lift distributions
xDelta = 100e-3;
fuselageWidth = 2;
x = 0:xDelta:b/2;
cy = cr - 2/b*(cr-ct)*x; % trapezoid. And local chord

   
%% volume in wing - NASA1015

% start at (0,0)
FX_74_Cl5_140 = [1.000000  0.000090
  0.997540  0.001890
  0.990700  0.006240
  0.980370  0.011500
  0.966980  0.017410
  0.950440  0.024020
  0.930640  0.031880
  0.907750  0.040010
  0.882020  0.048480
  0.853700  0.057180
  0.823090  0.066160
  0.790480  0.075490
  0.756160  0.085100
  0.720430  0.094820
  0.683590  0.104460
  0.645940  0.113850
  0.607780  0.122870
  0.569370  0.131370
  0.530990  0.139160
  0.492650  0.146040
  0.454350  0.151770
  0.416380  0.156060
  0.378870  0.158680
  0.342040  0.159440
  0.306090  0.158200
  0.271200  0.154930
  0.237600  0.149640
  0.205490  0.142430
  0.175040  0.133440
  0.146480  0.122920
  0.119990  0.111100
  0.095760  0.098260
  0.073950  0.084590
  0.054680  0.070300
  0.038110  0.055760
  0.024330  0.041450
  0.013380  0.027690
  0.005480  0.015180
  0.000980  0.005180
  0.000000 -0.000210
  0.000980 -0.004350
  0.005480 -0.007870
  0.013380 -0.008710
  0.024330 -0.007540
  0.038110 -0.005390
  0.054680 -0.002920
  0.073950 -0.000220
  0.095760  0.002700
  0.119990  0.005840
  0.146480  0.009210
  0.175040  0.012790
  0.205490  0.016510
  0.237600  0.020300
  0.271200  0.024100
  0.306090  0.027860
  0.342040  0.031520
  0.378870  0.035030
  0.416380  0.038320
  0.454350  0.041340
  0.492650  0.044000
  0.530990  0.046240
  0.569370  0.048010
  0.607780  0.049250
  0.645940  0.049920
  0.683590  0.049970
  0.720430  0.049360
  0.756160  0.048070
  0.790480  0.046080
  0.823090  0.043340
  0.853700  0.039780
  0.882020  0.035370
  0.907750  0.030130
  0.930640  0.024340
  0.950440  0.018510
  0.966980  0.012900
  0.980370  0.007560
  0.990700  0.003190
  0.997540  0.000680
  1.000000 -0.000030];

FX_72_MS_150A = [1.000000  0.000000
  0.998930  0.000370
  0.990390  0.003310
  0.973470  0.009140
  0.948440  0.017750
  0.915730  0.029010
  0.875920  0.042720
  0.853550  0.050420
  0.829670  0.058640
  0.804380  0.067350
  0.777790  0.076510
  0.750000  0.086070
  0.721140  0.095170
  0.691340  0.104250
  0.660720  0.112940
  0.629410  0.121370
  0.597550  0.129160
  0.565260  0.136310
  0.532700  0.142470
  0.500000  0.147810
  0.467300  0.151750
  0.434740  0.154440
  0.402450  0.155100
  0.370590  0.154340
  0.339280  0.151710
  0.308660  0.148130
  0.278860  0.143250
  0.250000  0.137530
  0.222210  0.130840
  0.195620  0.123390
  0.170330  0.115280
  0.146450  0.106570
  0.124080  0.097330
  0.103320  0.087730
  0.084270  0.077850
  0.066990  0.067770
  0.051560  0.057670
  0.038060  0.047690
  0.026530  0.037910
  0.017040  0.028670
  0.009610  0.019850
  0.004280  0.012520
  0.001070  0.006790
  0.000000  0.000000
  0.001070 -0.005960
  0.004280 -0.010560
  0.009610 -0.012460
  0.017040 -0.013800
  0.026530 -0.014320
  0.038060 -0.014860
  0.051560 -0.014850
  0.066990 -0.014610
  0.084270 -0.014070
  0.103320 -0.013280
  0.124080 -0.012310
  0.146450 -0.011130
  0.170330 -0.009750
  0.195620 -0.008260
  0.222210 -0.006590
  0.250000 -0.004800
  0.278860 -0.002770
  0.308660 -0.000720
  0.339280  0.001710
  0.370590  0.004230
  0.402450  0.007410
  0.434740  0.010900
  0.467300  0.015060
  0.500000  0.018810
  0.532700  0.022660
  0.565260  0.026090
  0.597550  0.029240
  0.629410  0.031700
  0.660720  0.033630
  0.691340  0.034730
  0.721140  0.035220
  0.750000  0.034870
  0.777790  0.033910
  0.804380  0.032230
  0.829670  0.030050
  0.853550  0.027350
  0.875920  0.024370
  0.915730  0.017930
  0.948440  0.011600
  0.973470  0.006170
  0.990390  0.002360
  0.998930  0.000340
  1.000000  0.000000];

zeroIndex = find(FX_72_MS_150A(:,1) == 0);
aerofoilX      = FX_72_MS_150A(zeroIndex:end,1);
aerofoilTop    = flip(FX_72_MS_150A(1:zeroIndex,2));
aerofoilBottom = FX_72_MS_150A(zeroIndex:end,2);


figure(1)
plot(aerofoilX,aerofoilTop,'b', aerofoilX,aerofoilBottom,'b')
    

% volume required
% for PL = 2500kg, endurance = 50hr
f_used = sum(fused); % kg
rho_fuel = 804; % kg/m^3
         
% sums area inside aerofoil, then scales x and y by multiplying by local
% chord c, and gets volume of element with xDelta  
wingArea = sum((aerofoilTop(1:end-1)-aerofoilBottom(1:end-1)).*diff(aerofoilX));

wingFuelFactorY = zeros(1,length(x));
frontSparLoc = 0.1;
rearSparLoc = 0.65;
foldLoc = 5.2/(b/2);
wingFuelFactorY(round(length(x)*frontSparLoc) : round(length(x)*rearSparLoc)) = 1;
wingFuelFrac = sum(wingFuelFactorY)/length(wingFuelFactorY);
fuelStopIndex = round(length(x)*foldLoc*0.9); % foldLoc% of length used for fuel (assume not all of this volume for fuel (eg. ribs) so take 90%) (array index for last part of wing length used for fuel)
wingFuelFactorX = wingFuelFrac * [ones(1,fuelStopIndex), zeros(1,length(x)-fuelStopIndex)]; % 70% of fuel useable in wing area (

wingFuelWeight = wingArea*cy.^2*rho_fuel*xDelta .* wingFuelFactorX; % multiply by useable space in fuel % c^2 to scale x and y directions to chord.
volumeWing = sum(wingFuelWeight)/rho_fuel;

fprintf('Volume in 1 wing: %.3f m^3\n', volumeWing);

V_required = f_used/rho_fuel;
fprintf('Volume required: %.3f m^3\n', V_required);
internalFuel = V_required - 2*volumeWing;
internalFuelWeight = internalFuel*rho_fuel; 
fprintf('Fuel Amount left (wings %.0f%% full): %.3f m^3\n', wingFuelFrac*100, internalFuel);

%% Material list
% AL 7075-T6 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T6
Materials.AL7075.E = 71.7e9; % Pa 
Materials.AL7075.rho = 2810; % (kg/m^3)
Materials.AL7075.shearStrength = 331e6; % Pa
Materials.AL7075.shearModulus = 26.9e9; % Pa
Materials.AL7075.ultimateTensile = 572e6; % Pa
Materials.AL7075.tensileYield = 503e6; % Pa

% AL 7075-T73 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T73
Materials.AL7075T73.E = 72e9; % Pa 
Materials.AL7075T73.rho = 2810; % (kg/m^3)
Materials.AL7075T73.shearStrength = 300e6; % Pa
Materials.AL7075T73.shearModulus = 26.9e9; % Pa
Materials.AL7075T73.ultimateTensile = 505e6; % Pa
Materials.AL7075T73.tensileYield = 435e6; % Pa

% AL 7475-T61 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7475T61
Materials.AL7475.E = 70.3e9; % Pa 
Materials.AL7475.rho = 2810; % (kg/m^3)
Materials.AL7475.shearStrength = 330e6; % Pa
Materials.AL7475.shearModulus = 27e9; % Pa
Materials.AL7475.ultimateTensile = 565e6; % Pa
Materials.AL7475.tensileYield = 490e6; % Pa

% AL 7475-T61 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7475T61
Materials.AL7475.E = 70.3e9; % Pa 
Materials.AL7475.rho = 2810; % (kg/m^3)
Materials.AL7475.shearStrength = 330e6; % Pa
Materials.AL7475.shearModulus = 27e9; % Pa
Materials.AL7475.ultimateTensile = 565e6; % Pa
Materials.AL7475.tensileYield = 490e6; % Pa

% AL 7475-T7651 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7475T765
Materials.AL7475T7651.E = 70.3e9; % Pa 
Materials.AL7475T7651.rho = 2810; % (kg/m^3)
Materials.AL7475T7651.shearStrength = 330e6; % Pa
Materials.AL7475T7651.shearModulus = 27e9; % Pa
Materials.AL7475T7651.ultimateTensile = 565e6; % Pa
Materials.AL7475T7651.tensileYield = 490e6; % Pa

% AL2024 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T4
Materials.AL2024.E = 73.1e9; % Pa 
Materials.AL2024.rho = 2780; % (kg/m^3)
Materials.AL2024.shearStrength = 283e6; % Pa
Materials.AL2024.shearModulus = 28e9; % Pa
Materials.AL2024.ultimateTensile = 469e6; % Pa
Materials.AL2024.tensileYield = 324e6; % Pa

% Titanium Ti-6Al-4V (Grade 5), Annealed - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MTP641
Materials.TI64.E = 113.8e9; % Pa 
Materials.TI64.rho = 4430; % (kg/m^3)
Materials.TI64.shearStrength = 550e6; % Pa
Materials.TI64.shearModulus = 44e9; % Pa
Materials.TI64.ultimateTensile = 950e6; % Pa
Materials.TI64.tensileYield = 880e6; % Pa


%% force, moment calculation

% adjust for other flight conditions with differnt MTOW and load factor
% Cant do N=2.5 at MTOW but can do N=1But wont be doign 2.5 at MTOW


N = 2.5; % load factor
W = N*22500; % (kg)
rho = 1.225; %*3.468/14.696; % factor for 35000ft
V = sqrt(N*TOW*g/(0.5*rho*clmax*S)); % knots

Ltrap = 2*W/(b*(1+taper_r))*(1-2*x/b*(1-taper_r));
Lellip = 4*W/(pi*b)*sqrt(1-(2*x/b).^2);
Lavg = (Ltrap + Lellip)/2; % kg/m
LPerPanel = Lavg*xDelta; % kg


M0 = 120/340; % max flight at sea level
tc_ratio = 0.15;
wingWeightCalc = 0.453592 * 1.75 * 0.00428*(S*0.3048)^0.48*(AR*M0^0.43*(W*g*4.44822*N)^0.84*taper_r^0.14)/((100*tc_ratio)^0.76*(cosd(w_sweep))^1.54); % (N) - Subsonic Aircraft - Nicolai
wWPerSpan = wingWeightCalc/S * cy; % (kg/m) - wing weight per span
wWPerSpan(round(foldLoc*length(x))) = wWPerSpan(round(foldLoc*length(x))) + 500; % extra weigth at fold location for everything there


payloadLocation = [0.1 foldLoc]; % percent of wing % this number of payload on each wing
payload = 0;%3500/2; % worst case is no payload because then greater net vertical force
payloadPerSpan = zeros(1,length(x));
for i = 1:length(payloadLocation)
    payloadIndex = round(payloadLocation(i)*length(x));
    devP = round(0.1/xDelta); % payload from -devPayload to +devPayload around point load (0.1m either side)
    bound = -devP+payloadIndex:devP+payloadIndex;
    payloadPerSpan(bound) = payload/length(bound); % kg/m  
%     payloadWeight(payloadIndex) = payload/length(payloadLocation); % kg/m  (point load) --------- should I size it for 2 payloads or one on outer edge? i.e. two half loads or 1 full load
end
payloadY = 0.2; % m - distance from wing to payload
payloadDrag = Cdo_pod * 0.5 * rho * V^2 * S / g; % kg - Drag on 1 payload

fuelPerSpan = wingArea*cy.^2*rho_fuel .* wingFuelFactorX; % kg/m
selfWeight = N*(wWPerSpan + payloadPerSpan + fuelPerSpan); % kg/m


nYDelta = 100;
for i = 1:length(x)
    c_mat(:,i) = linspace(0,cy(i),nYDelta);
%     q_mat(:,i) = 2 * Lavg(i)/c(i) * (1 - c_mat(1:end-1,i)/c(i)); % kg/m^2 --- triangular chordwise distribution
    
    x015 = 0.15*cy(i); % chordwise length at 0.15 chord
    [~, xi015] = min(abs(x015 - c_mat(:,i))); % index closest to 0.15c position
    cModified = [0.15*cy(i)*ones(xi015,1); c_mat(xi015+1:end,i)]; % modidy leading edge values for next equation. Constant value at leading edge to produce flat part
    q_mat(:,i) = 40/23*Lavg(i)/cy(i) * (1 - 1/0.85 * (cModified/cy(i) - 0.15)); % triangular chordwise distribution with flat leading edge as per Raymer pg 345
    
end
q_mat(:,end) = 0;

%q = 0.5*c.*q; % to account for triangular chordwise lift distribution - combine total force at spanwise position into 1 value to use as N/m even though units would be N
selfWeightMat = repmat(selfWeight./cy,size(q_mat,1),1); % kg/m^2
q = q_mat - selfWeightMat; % kg/m^2
q_orig = Lavg; % kg/m
q_mat_orig = q_mat; % kg/m^2
%totalLift = sum(q_orig*xDelta); % assume triangular shaped distribution along chord (largest load at leading edge)
totalLift = sum(sum(q_mat.*cy/nYDelta)*xDelta);
fprintf('total Lift (half): %.0f kg\n', totalLift);
fprintf('weight (half): %.0f kg\n\n', W/2);

LPerSpanOrig = sum(q_mat.*cy/nYDelta);
LPerSpan  = sum(q.*cy/nYDelta);
LPerChord = sum(q'*xDelta);
figure(6)
plot(LPerChord);


h = 0.15*cy; % m - height of beam at each section
b_cap0 = 111e-3; % m - constant beam width (value for plot)
t_cap0 = 10e-3; % m 

% % FX_74_Cl5_140
% b_cap1 = 142e-3; % m
% t_cap1 = 10e-3; % m
% b_cap2 = 111e-3; % m
% t_cap2 = 10e-3; % m

% FX_72_MS_150A
b_cap1 = 230e-3; % m - front spar
t_cap1 = 10e-3; % m
b_cap2 = 180e-3; % m - rear spar % need deflection of both beams to be the same so there is no twisting from normal lift force
t_cap2 = 10e-3; % m

b_cap_vec = b_cap2; % m
t_cap_vec = t_cap2; % m

FS = 1.5; % STANAG Subpart C 303 - "F.O.S no lower than 1.5 for structures whose failure would lead to a hazardous or more serious failure condition"

[~, frontIndex] = min(abs(frontSparLoc - aerofoilX)); % index closest to front spar location
[~, rearIndex] = min(abs(rearSparLoc - aerofoilX)); % index closest to rear spar location
h1 = (aerofoilTop(frontIndex)-aerofoilBottom(frontIndex)) * cy;
h2 = (aerofoilTop(rearIndex)-aerofoilBottom(rearIndex)) * cy;

% assign material values
% beamUsed = 'AL7075';

disp('Weight of half beam')
materials = fieldnames(Materials);
for mIndex = 4%:length(materials)
        
    beamUsed = materials{mIndex};
    E = Materials.(beamUsed).E;
    rhoBeam = Materials.(beamUsed).rho; % (kg/m^3)
    shearStrength = Materials.(beamUsed).shearStrength; % Pa
    shearModulus = Materials.(beamUsed).shearModulus; % Pa
    ultimateTensile = Materials.(beamUsed).ultimateTensile; % Pa
    tensileYield = Materials.(beamUsed).tensileYield; % Pa
    
    m = 1;
    n = 1;
    
for m = 1:length(t_cap_vec)
    t_cap2 = t_cap_vec(m);
    
    for n = 1:length(b_cap_vec)
        b_cap2 = b_cap_vec(n);
        
        beamWeightEq = @(b_cap,t_cap,h) rhoBeam*t_cap*(2*b_cap+h-2*t_cap); % kg/m
        beamWeight = beamWeightEq(b_cap1,t_cap1,h1) + beamWeightEq(b_cap2,t_cap2,h2);

        %{
        %bending
        https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10.pdf
        
        %force distribution: front - rear spar (for high aspect ratio wing)
        https://pdfs.semanticscholar.org/e288/128ae125e0a6b8dd949d7e1afa1b97247bb0.pdf
         
        %beam stiffness in bending: k = 3EI/L^3
        https://engineering.sjsu.edu/e10/wp-content/uploads/Structure_Stiffness_S13.pdf
        %} 
        
        % moment of inertia
        I_cap = @(b_cap,t_cap,h) b_cap*t_cap.^3/12 + b_cap*t_cap*(h/2-t_cap/2).^2;
        I_middle = @(b_cap,t_cap,h) b_cap*(h-2*t_cap).^3/12;
        I_combined = @(b_cap,t_cap,h) 2*I_cap(b_cap,t_cap,h) + I_middle(b_cap,t_cap,h);
        
        A_beamEq = @(b_cap,t_cap,h) 2*b_cap*t_cap + b_cap*(h-2*t_cap);
        
        I1 = I_combined(b_cap1,t_cap1,h1);
        I2 = I_combined(b_cap2,t_cap2,h2);
        A_beam1 = A_beamEq(b_cap1,t_cap1,h1);
        A_beam2 = A_beamEq(b_cap2,t_cap2,h2);
        K1 = 3*E*I1/(b/2)^3; % front spar stiffness
        K2 = 3*E*I2/(b/2)^3; % rear spar stiffness
%         P1 = q .* K1./(K1+K2);
%         P2 = q .* K2./(K1+K2);
%         P1(end) = 0;  
%         P2(end) = 0;
        P1 = q .* (rearSparLoc-0.25)/(rearSparLoc-frontSparLoc); % if distance determines force distribution
        P2 = q .* (0.25-frontSparLoc)/(rearSparLoc-frontSparLoc);
        
        
% spanwise forces/moments
        figure(2)
        plot(x,sum(q_mat.*cy/nYDelta), x,sum(selfWeightMat.*cy/nYDelta), x,sum(q_mat.*cy/nYDelta)-sum(selfWeightMat.*cy/nYDelta), x,sum(P1.*cy/nYDelta), x,sum(P2.*cy/nYDelta))
        ylabel('Force/span (kg/m)')
        xlabel('span location (m)')
        legend('Lift','Self Wing Weight','Net Vertical', 'Front','Rear')
        [Sf1, Sf2]       = deal(zeros(1,length(x))); % kg - shear force
        [M1, M2]         = deal(zeros(1,length(x))); % kgm - bending moment
        [Ss1, Ss2]       = deal(zeros(1,length(x))); % kg/m^2 - shear stress
        [theta1, theta2] = deal(zeros(1,length(x))); % rad - deflection angle
        [w1, w2]         = deal(zeros(1,length(x))); % m - deflection

        [SsF1,SsF2,SsM1,SsM2] = deal(zeros(1,length(x)));
        for i = length(x)-1:-1:1 % 0 stress and moment at wing tip
            Sf1(i) = Sf1(i+1) - 0.5*(sum(P1(:,i+1))+sum(P1(:,i)))*xDelta*(cy(i)/nYDelta); % kg
            Sf2(i) = Sf2(i+1) - 0.5*(sum(P2(:,i+1))+sum(P2(:,i)))*xDelta*(cy(i)/nYDelta); % kg
            
            M1(i) = M1(i+1) - 0.5*(Sf1(i+1)+Sf1(i))*xDelta; % kgm
            M2(i) = M2(i+1) - 0.5*(Sf2(i+1)+Sf2(i))*xDelta;
            
            SsF1(i) = Sf1(i)./A_beam1(i); % transverse stress
            SsM1(i) = M1(i).*(h1(i)/2)./I1(i); % bending stress
            SsF2(i) = Sf2(i)./A_beam2(i);
            SsM2(i) = M2(i).*(h2(i)/2)./I2(i); 
            Ss1(i)  = abs(SsF1(i) - SsM1(i)); % total shear stress
            Ss2(i)  = abs(SsF2(i) - SsM2(i));
        end
        for j = 2:length(x) % 0 deflection at wing root
            theta1(j) = theta1(j-1) + 0.5*(M1(j)/E/I1(j) + M1(j-1)/E/I1(j-1))*xDelta; % non-dimensional 
            theta2(j) = theta2(j-1) + 0.5*(M2(j)/E/I2(j) + M2(j-1)/E/I2(j-1))*xDelta;
            w1(j) = w1(j-1) + 0.5*(theta1(j)+theta1(j-1))*xDelta; % m 
            w2(j) = w2(j-1) + 0.5*(theta2(j)+theta2(j-1))*xDelta;
        end

        
        
% chordwise forces and moments
        Mpayload = payloadY*payloadDrag * length(payloadLocation);

        
%         if t_cap == t_cap2 && b_cap == b_cap2
            figure(mIndex + 2)
            subplot(2,2,1)
            %subplot(4,1,1)
            plot(x,Sf1, x,Sf2) %, 'DisplayName', beamUsed)
            ylabel('Shear Force (N)')
            xlabel('span location (m)')
            %title(['b_{cap} = ',num2str(b_cap*1000), ' mm, t_{cap} = ',num2str(t_cap*1000),' mm'])
            ax = gca;
            ax.TitleFontSizeMultiplier = 0.6;
            %legend_labels{mIndex} = beamUsed;
            title(beamUsed);
            xlim([0 b/2])
            hold on
            
            subplot(2,2,2)
            %subplot(4,1,3)
            plot(x,Ss1*FS, x,Ss2*FS, [0 b/2],[shearStrength shearStrength]/g,'--g')  
            ylabel('Shear Stress (N/m^2)')
            xlabel('span location (m)')
            legend('front', 'rear', 'max allowed with F.S', 'location','SE')
            xlim([0 b/2])
            hold on
            
            subplot(2,2,3)
            %subplot(4,1,2)
            plot(x,M1, x,M2)
            ylabel('Bending Moment (Nm)')
            xlabel('span location (m)')
            xlim([0 b/2])
            hold on
            
            subplot(2,2,4)
            plot(x,rad2deg(theta1), x,rad2deg(theta2))
            ylabel('Angular Deflection (deg)')
            xlabel('span location (m)')
            xlim([0 b/2])
            hold on
            
%             subplot(2,2,4)
%             %subplot(4,1,4) %if deflection of both beams is the same, then no twisting from lift distribution
%             plot(x,w1, x,w2)
%             ylabel('Deflection (m)')
%             xlabel('span location (m)')
%             xlim([0 b/2])
%             hold on
%         end
        
        % margin of safety
        MS1(m,n) = (shearStrength/g)/(max(Ss1)*FS) - 1; % divide by 'g' since Sf in 'kg'
        MS2(m,n) = (shearStrength/g)/(max(Ss2)*FS) - 1;
        
        maxS(m,n) = min(Sf1);
        maxM(m,n) = max(M1);
        maxW(m,n) = max(w1);
        
        fprintf('%12s weight: %.0f kg | %.0f, %.0f, %.0f, %.0f | %.4f, %.4f\n', beamUsed,sum(beamWeight*xDelta),b_cap1*1000,t_cap1*1000,b_cap2*1000,t_cap2*1000,MS1(m,n),MS2(m,n));
    end
end
end


%legend('show');
%legend(legend_labels,'location','NW')

% for plotting variations in I beam dimensions
% figure(5)
% subplot(2,1,1)
% for j = 1:size(maxS,1)
%     plot(b_cap_vec*1000, MS1(j,:),'linewidth',3) %maxS(j,:),'linewidth',3)
%     legend_labels{j} = sprintf('t_{cap} = %.0f mm.', t_cap_vec(j)*1000);
%     hold on;
% end
% legend(legend_labels)
% xlabel('b cap width (mm)');
% %ylabel('min shear (N/m^2)');
% ylabel('Margin Of Safety')
% 
% subplot(2,1,2)
% for j = 1:size(maxM,1)
%     plot(b_cap_vec*1000,MS2(j,:),'linewidth',3) %maxM(j,:),'linewidth',3)
%     %legend_labels{j} = sprintf('t_{cap} = %.2f mm.', t_cap_vec(j));
%     hold on;
% end
% legend(legend_labels)
% xlabel('b cap width (mm)');
% %ylabel('max moment (Nm)');
% ylabel('Margin Of Safety')

% subplot(2,2,3)
% for j = 1:size(maxW,1)
%     plot(b_cap_vec*1000,maxW(j,:),'linewidth',3)
%     legend_labels{j} = sprintf('t_{cap} = %.2f m.', t_cap_vec(j));
%     hold on;
% end
% legend(legend_labels)
% xlabel('b cap width (mm)');
% ylabel('max deflection (m)');