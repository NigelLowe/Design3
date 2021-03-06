clear all, 
close all,
clc

% DESIGN 3 - Preliminary Sizing - Group 3

% Call Parameters
albatross_parameters_maritime;
%albatross_parameters_airfield;

%Mission Req
% Mission 1 MAX PL and 40 hrs En
mission_1_x = [40 40];
mission_1_y = [0 60e4];

% Mission 2 half PL and 50 hrs En
mission_2_x = [50 50];
mission_2_y = [0 60e4];

%Vectors of PL and Endurance
pl_vector = [500, 3500]; %linspace(0,3500,3);
en_vector = linspace(10,50,9); 


%Wait bar for sanity
hwait = waitbar(0,'Please wait...');
h_bar = 1;
for m_index = 1:length(pl_vector)
    
    PL = pl_vector(m_index);
    for n_index = 1:length(en_vector)
        
        endurance_target = en_vector(n_index);
        
        climb_cruise_iteration
        
        TOW_results(m_index,n_index) = TOW;

        % Wait bar
        hwait = waitbar(h_bar/(length(pl_vector) * length(en_vector)));
        h_bar = h_bar + 1;
    end
end

close(hwait)

%% Plot lines of contant Payload
clf
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
assumptions = sprintf('cdo: %.4f, TSFC: %.6f, e: %.2f, alt: %.0f ft, S: %.0f m^2, AR: %.0f, v cruise %.0f m/s, clmax: %0.f, P: %0.f kW', cdo,TSFC,e,cruise_alt*3.281,S,AR,v_cruise,clmax,P/1000);
annotation('textbox',...
    [0.55 0.86 0.35 0.05],...
    'String',assumptions,...
    'LineStyle','none',...
    'FontSize',12,...
    'FitBoxToText','off');

tow = TOW_results(~isinf(TOW_results));
ylim([max(min(min(tow))-1e4,0),max(max(tow))+1e4])
xlim([10,50])

legend(legend_labels,'location','NW')


%% Take off distance

for m_index = 1:length(pl_vector)
    for n_index = 1:length(en_vector)
        
        m = TOW_results(m_index,n_index);
        Vstall = sqrt(m*g/(0.5*clmax*rho0*S));
        Vr = 1.1*Vstall;
        V = Vr/sqrt(2);

        T = P*prop_n/V;

        cd = cd0 + k*(clmax-clmin)^2;
        D = 0.5*cd*rho0*V^2*S;
        L = 0.5*clmax*rho0*V^2*S;
        
        mu = 0.02;
        F = mu*(m*g - L);

        a = (T - D - F)/m;

        s_dist(m_index,n_index) = Vr^2/(2*a);

    end
end

figure(2);
hold on
grid on
ylabel('Take-Off Distance (m.)')
xlabel('Endurance (hrs.)')
set(gca,'FontSize',18)

% Plot lines of constant PL
for j = 1:size(s_dist,1)
    
    plot(en_vector,s_dist(j,:),'linewidth',3)
    legend_labels{j} = sprintf('%.0f kg.', pl_vector(j));
    
end

%assumptions box
assumptions = sprintf('cdo: %.4f, TSFC: %.6f, e: %.2f, alt: %.0f ft, S: %.0f m^2, AR: %.0f', cdo,TSFC,e,cruise_alt*3.281,S,AR);
annotation('textbox',...
    [0.33 0.86 0.21 0.05],...
    'String',assumptions,...
    'LineStyle','none',...
    'FontSize',12,...
    'FitBoxToText','off');

ylim([0,min(1000,max(max(s_dist))*1.05)])
xlim([10,50])

legend(legend_labels,'location','NW')
