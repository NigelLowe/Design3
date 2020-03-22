clear all, 
close all,
clc

% DESIGN 3 - Preliminary Sizing - Group 3

%Wing + Aero Parameters
S   = 135;        %m2 % -------------------
AR  = 40;             % -------------------
b   = sqrt(S*AR); %m
c   = S/b;

e   = 0.7;
cdo = 0.03;
k   = 1/(pi*e*AR);

%Engine
TSFC = 14*1e-6;   %kg/(s.N) % -------------------

%mission parameters
reach_toc        = 6;  %hrs
cruise_alt       = 40000/3.281; %m
v_cruise         = 180*0.51444; %m/s

target_roc       = cruise_alt/reach_toc/3600; %m/s
% should replace this with specific excess power
% roc = arcsin( (dh/dt) / V );
ld_climb         = 1/(2*sqrt(k*cdo)); %min l_d
cl_climb         = sqrt(pi*AR*cdo*e);
% is this for climb or decent?

reserves         = 0.01; %decent/approach/taxi % -------------------

clmax = 1.8;
clmin = 0.2; % -------------------


%Mission Req
% Mission 1 MAX PL and 40 hrs En
mission_1_x = [40 40];
mission_1_y = [0 60e4];

% Mission 2 half PL and 50 hrs En
mission_2_x = [50 50];
mission_2_y = [0 60e4];

%Vectors of PL and Endurance
pl_vector = linspace(0,3500,3);
en_vector = linspace(40,60,10); 


%Wait bar for sanity
hwait = waitbar(0,'Please wait...');
h_bar = 1;
for m_index = 1:length(pl_vector)
    
    PL = pl_vector(m_index);
    for n_index = 1:length(en_vector)
        
        endurance_target = en_vector(n_index);
        
        climb_cruise_iteration
        
        TOW_results(m_index,n_index) = TOW;
        Power_required(m_index,n_index) = Power;
        
        % Wait bar
        hwait = waitbar(h_bar/(length(pl_vector) * length(en_vector)));
        h_bar = h_bar + 1;
    end
end

close(hwait)

%% Plot lines of contant Payload
fig = figure(1);
hold on
grid on
ylabel('TOW (kg.)')
xlabel('Endurance (hrs.)')
set(gca,'FontSize',18)

% Plot lines of constant PL
for j = 1:size(TOW_results,1)
    
    plot(en_vector,TOW_results(j,:),'linewidth',3)
    legend_labels{j} = sprintf('%.0f kg.', pl_vector(j));
    
end

%assumptions box
assumptions = sprintf('cdo: %.4f, TSFC: %.6f, e: %.2f, alt: %.0f ft, S: %.0f m^2, AR: %.0f', cdo,TSFC,e,cruise_alt*3.28084,S,AR);
annotation('textbox',...
    [0.334375 0.860851505711319 0.208854166666667 0.0492253374870202],...
    'String',assumptions,...
    'LineStyle','none',...
    'FontSize',12,...
    'FitBoxToText','off');

% Plot mission lines
plot(mission_1_x,mission_1_y,'--','color','k')
plot(mission_2_x,mission_2_y,'--','color','k')

ylim([1e4,max(max(TOW_results))+1e4])
xlim([39,61])

plot([39.5 60.5],[3.7e4 3.7e4],...
   'color','k','linewidth',3)

legend(legend_labels,'location','NW')

figure(2);
hold on
grid on
ylabel('Power (kW)')
xlabel('Endurance (hrs.)')
set(gca,'FontSize',18)

% Plot lines of constant PL
for j = 1:size(Power_required,1)
    
    plot(en_vector,Power_required(j,:)/1000,'linewidth',3)
    legend_labels{j} = sprintf('%.0f kg.', pl_vector(j));
    
end
%% Take off distance
% use weight for 50hrs
W = TOW_results(1,6);
g = 9.81;
rho = 1.225;
CL = 1.8; % CLmax
CLmin = 0.1;

Vstall = sqrt(W*g/(0.5*CL*rho*S));
Vr = 1.1*Vstall;
V = Vr/sqrt(2);

n = 1;
P = n*8202000; % 1 TP-400 engine
eff = 0.85;
T = P*eff/V;

Cd0 = 0.03;
Cd = Cd0 + k*(CL-CLmin)^2;
D = 0.5*Cd*rho*V^2*S;

L = 0.5*CL*rho*V^2*S;
mu = 0.02;
F = mu*(W*g - L);

a = (T - D - F)/W;

s = Vr^2/(2*a);
disp(['takeoff dist: ',num2str(round(s)),' m']);