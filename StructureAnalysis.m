% Structural Analysis
% 460368355

if ~exist('plotOtherGraphs','var') % if statement for this file use in other functions
clear
    
    % general parameters
    albatross_parameters_maritime;
    %albatross_parameters_airfield;
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
    pl_num = 3500;
    en_num = 24;
end
plotOtherGraphs = 'no';
prelim_report_code;
clc
close all

% flight condition
V = 80; % m/s

% for taper ratio input
taper_r = 0.45; % from graph in lecture 5 slide 19. Minimum drag at 0 sweep
% predator - taper_r = 0.25
% reaper - taper_r = 0.8
% global hawk - taper_r = 0.3
%{ 
https://www.semanticscholar.org/paper/A-Genetic-Algorithm-Incorporating-Design-Choice-for-Mull/0742bbd9eb7ea5651c01241ae2357bde480a4625?tab=abstract&citingPapersSort=is-influential&citingPapersLimit=10&citingPapersOffset=0&year%5B0%5D=&year%5B1%5D=&citedPapersSort=is-influential&citedPapersLimit=10&citedPapersOffset=10
%}

cr = 2*b/(AR*(1+taper_r)); % m - chord at tip
ct = cr * taper_r; % m - chord at root
mean_ac = 2/3 * cr * (1+taper_r+taper_r^2) / (1+taper_r);


% load at center (mean line method) for lift distribution
% mean line between ellipse and trapezoid with same area.
z = 2/pi * (ct + cr); % N
% ellipse equation: (x/(b/2))^2 + (y/z)^2 = 1;

% equations of lift distributions
xDelta = 0.1e-3;
x = 0:xDelta:b/2;
yt = cr - 2/b*(cr-ct)*x; % trapezoid. And local chord
ye = z*sqrt(1 - (2*x/b).^2); % ellipse
c = (yt + ye)/2; % mean
c(end) = 0;
% figure(1)
plot(x,yt, x,ye, x,c)
legend('trapezoid','ellipse','mean')
ylabel('lift distribution shape')
xlabel('span location (m)')
xlim([0 b/2])
   
%% volume in wing - NASA1015

% start at (0,0)
nasa1015x = [0.000000 ...
             0.001621 0.006475 0.014529 0.025732 0.040010 0.057272 0.077405 0.100279 0.125745 0.153638 ...
             0.183777 0.215968 0.250000 0.285654 0.322698 0.360891 0.399987 0.439732 0.479867 0.520133 ...
             0.560268 0.600013 0.639109 0.677302 0.714346 0.750000 0.784032 0.816223 0.846362 0.874255 ...
             0.899721 0.922595 0.942728 0.959990 0.974268 0.985471 0.993525 0.998379 1.000000];
         
nasa1015top = [0.000000 0.000577 0.002244 0.004834 0.008107 0.011820 0.015794 0.020019 0.024801 0.030724 ...
             0.037889 0.046102 0.055549 0.066230 0.077492 0.088455 0.098407 0.106947 0.113813 0.118919 ...
             0.122351 0.124264 0.124809 0.124114 0.122299 0.119458 0.115663 0.110974 0.105433 0.099088 ...
             0.092004 0.084202 0.075746 0.066677 0.057142 0.046990 0.036312 0.026197 0.017070 0.000000];
nasa1015top = flip(nasa1015top);

nasa1015bottom = [0.000000 ...
             -0.001976 -0.004930 -0.007513 -0.010010 -0.012525 -0.014983 -0.017226 -0.019318 -0.021214 ...
             -0.022877 -0.024304 -0.025466 -0.026357 -0.026977 -0.027302 -0.027330 -0.027065 -0.026505 ...
             -0.025652 -0.024494 -0.023028 -0.021274 -0.019239 -0.016865 -0.014081 -0.010938 -0.007663 ...
             -0.004646 -0.002130 -0.000215  0.001069  0.001761  0.001957  0.001792  0.001378  0.000884 ...
             0.000429 0.000113 0.000000];

% figure(2)
% plot(nasa1015x,nasa1015top,'b', nasa1015x,nasa1015bottom,'b')
    

% volume required
% for PL = 2500kg, endurance = 50hr
f_used = sum(fused); % kg
rho_fuel = 804; % kg/m^3
         
% sums area inside aerofoil, then scales x and y by multiplying by local
% chord c, and gets volume of element with xDelta  
wingArea = sum((nasa1015top(1:end-1)-nasa1015bottom(1:end-1)).*diff(nasa1015x));
%volumeWing = sum(wingArea * c.^2*xDelta) * 0.7*0.4; % 70% of wing available for fuel usage, 40% length of wing used (no fuel in folding part)
%wingFuelWeight = volumeWing*rho_fuel; %

fuelStopIndex = round(length(x)*0.4); % 40% of length used for fuel (array index for last part of wing length used for fuel)
wingFuelFactor = 0.7 * [ones(1,fuelStopIndex), zeros(1,length(x)-fuelStopIndex)]; % 70% of fuel useable in wing area (
wingFuelWeight = wingArea*c.^2*rho_fuel*xDelta .* wingFuelFactor; % multiply by useable space in fuel
volumeWing = sum(wingFuelWeight)/rho_fuel;

fprintf('Volume in 1 wing: %.3f m^3\n', volumeWing);

V_required = f_used/rho_fuel;
fprintf('Volume required: %.3f m^3\n', V_required);
internalFuel = V_required - 2*volumeWing;
internalFuelWeight = internalFuel*rho_fuel; 
fprintf('Fuel Amount left (wings 70%% full): %.3f m^3\n', internalFuel);

%% Material list
% AL 7075-T6 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA7075T6
Materials.AL7075.E = 71.7e9; % Pa 
Materials.AL7075.rho = 2810; % (kg/m^3)
Materials.AL7075.shearStrength = 331e6; % Pa
Materials.AL7075.shearModulus = 26.9e9; % Pa
Materials.AL7075.ultimateTensile = 572e6; % Pa
Materials.AL7075.tensileYield = 503e6; % Pa

% AL2024 - http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=MA2024T4
Materials.AL2024.E = 73.1e9; % Pa 
Materials.AL2024.rho = 2780; % (kg/m^3)
Materials.AL2024.shearStrength = 283e6; % Pa
Materials.AL2024.shearModulus = 28e9; % Pa
Materials.AL2024.ultimateTensile = 469e6; % Pa
Materials.AL2024.tensileYield = 324e6; % Pa

%% force, moment calculation

n = 1; % g loading - highest net force at n <= 1
L = n * 0.5*rho0*V^2*clmax * c; % max load - sea level and max CL %%%%%%%%%%%%% dont involve weight in ligt calculation. I think it should
L_orig = L;
totalLift = 0.5 * sum(c.*L) * xDelta; % assume triangular shaped distribution along chord (largest load at leading edge)
fprintf('center Lift: %.0f N\n', L(1));
fprintf('total Lift: %.0f N\n', totalLift);
      
q = n*TOW*g/S * c; % way too low %%%%%%%% need to fix up
totalQ = 0.5 * sum(c.*q) * xDelta;
fprintf('total Lift Equation: %.0f N\n\n', totalQ);

h = 0.15*c; % m - height of beam at each section
b_cap0 = 0.15; % m - constant beam width (value for plot)
t_cap0 = 0.07; % m 

b_cap_vec = b_cap0; %0.05:0.01:0.3; % m
t_cap_vec = t_cap0; %0.01:0.02:0.1; % m

% assign material values
% beamUsed = 'AL7075';

materials = fieldnames(Materials);
for mIndex = 1:length(materials)
        
    beamUsed = materials{mIndex};
    E = Materials.(beamUsed).E;
    rhoBeam = Materials.(beamUsed).rho; % (kg/m^3)
    shearStrength = Materials.(beamUsed).shearStrength; % Pa
    shearModulus = Materials.(beamUsed).shearModulus; % Pa
    ultimateTensile = Materials.(beamUsed).ultimateTensile; % Pa
    tensileYield = Materials.(beamUsed).tensileYield; % Pa
    
    m = 1;
    n = 1;
    
% for m = 1:length(t_cap_vec)
    t_cap = t_cap_vec(m);
%     for n = 1:length(b_cap_vec)
    b_cap = b_cap_vec(n);
        
        beamWeight = rhoBeam*t_cap*(2*b_cap+h-2*t_cap)*xDelta;
        wingWeight = (wingFuelWeight + beamWeight)*g;
        L = L - wingWeight;
        figure(3)
        plot(x,L_orig, x,wingWeight, x,L, x,q)
        ylabel('Force (N)')
        xlabel('span location (m)')
        legend('Lift','Fuel Weight','Net Vertical','Equation')


        % bending
        %{
        https://ocw.mit.edu/courses/aeronautics-and-astronautics/16-01-unified-engineering-i-ii-iii-iv-fall-2005-spring-2006/systems-labs-06/spl10.pdf
        %}      
        Sh = zeros(1,length(x)); % N/m^2? - shear
        M = zeros(1,length(x)); % Nm - bending moment
        theta = zeros(1,length(x)); % rad - deflection angle
        w = zeros(1,length(x)); % m - deflection
        
        % moment of inertia
        I_cap = b_cap*t_cap.^3/12 + b_cap*t_cap*(h/2-t_cap/2).^2; 
        I_middle = b_cap*(h-2*t_cap).^3/12;
        I = 2*I_cap + I_middle;

        for i = length(x)-1:-1:1 % 0 stress and moment at wing tip
            Sh(i) = Sh(i+1) - 0.5*(0.5*L(i+1)*c(i+1)+0.5*L(i)*c(i))*xDelta;
            M(i) = M(i+1) - 0.5*(Sh(i+1)+Sh(i))*xDelta;
        end
        for j = 2:length(x) % 0 deflection at wing root
            theta(j) = theta(j-1) + 0.5*(M(j)/E/I(j) + M(j-1)/E/I(j-1))*xDelta;
            w(j) = w(j-1) + 0.5*(theta(j)+theta(j-1))*xDelta;
        end

        if t_cap == t_cap0 && b_cap == b_cap0
            figure(4)
            subplot(2,2,1)
            plot(x,Sh)%, 'DisplayName', beamUsed)
            ylabel('Shear (N)')
            xlabel('span location (m)')
            title(['b_{cap} = ',num2str(b_cap), ' m, t_{cap} = ',num2str(t_cap),' m'])
            ax = gca;
            ax.TitleFontSizeMultiplier = 0.8;
            xlim([0 b/2])
            hold on
            
            subplot(2,2,2)
            plot(x,M)
            ylabel('Bending Moment (Nm)')
            xlabel('span location (m)')
            xlim([0 b/2])
            hold on
            
            subplot(2,2,3)
            plot(x,rad2deg(theta))
            ylabel('Angular Deflection (ded)')
            xlabel('span location (m)')
            xlim([0 b/2])
            hold on
            
            subplot(2,2,4)
            plot(x,w)
            ylabel('Deflection (m)')
            xlabel('span location (m)')
            xlim([0 b/2])
            legend_labels{mIndex} = beamUsed;
            hold on
        end
        
        maxS(m,n) = min(Sh);
        maxM(m,n) = max(M);
        maxW(m,n) = max(w);
        
        fprintf('%s weight: %.0f kg\n', beamUsed,sum(beamWeight));
    %end
end

%legend('show');
legend(legend_labels,'location','NW')

% figure(5)
% subplot(2,2,1)
% for j = 1:size(maxS,1)
%     plot(b_cap_vec*1000,maxS(j,:),'linewidth',3)
%     legend_labels{j} = sprintf('t_{cap} = %.2f m.', t_cap_vec(j));
%     hold on;
% end
% legend(legend_labels)
% xlabel('b cap width (mm)');
% ylabel('min shear (N/m^2)');
% 
% subplot(2,2,2)
% for j = 1:size(maxM,1)
%     plot(b_cap_vec*1000,maxM(j,:),'linewidth',3)
%     legend_labels{j} = sprintf('t_{cap} = %.2f m.', t_cap_vec(j));
%     hold on;
% end
% legend(legend_labels)
% xlabel('b cap width (mm)');
% ylabel('max moment (Nm)');
% 
% subplot(2,2,3)
% for j = 1:size(maxW,1)
%     plot(b_cap_vec*1000,maxW(j,:),'linewidth',3)
%     legend_labels{j} = sprintf('t_{cap} = %.2f m.', t_cap_vec(j));
%     hold on;
% end
% legend(legend_labels)
% xlabel('b cap width (mm)');
% ylabel('max deflection (m)');