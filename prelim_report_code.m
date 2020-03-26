clear all, 

close all,

clc



% DESIGN 3 - Preliminary Sizing - Group 3



%Wing + Aero Parameters

S   = 68;        %m2

AR  = 18;            

b   = sqrt(S*AR); %m

c   = S/b;



e   = 0.85;

cdo = 0.03;

k   = 1/(pi*e*AR);



%Engine

TSFC = 13*1e-6;   %kg/(s.N)



%mission parameters

reach_toc        = 6;  %hrs

cruise_alt       = 40000/3.281; %m

loiter_point    = 4000e3; %4000 km

v_cruise        = convvel(250,'kts','m/s');
v_loiter        = convvel(200,'kts','m/s');

target_roc       = cruise_alt/reach_toc/3600; %m/s

ld_climb         = 1/(2*sqrt(k*cdo)); %min l_d

cl_climb         = sqrt(pi*AR*cdo*e);



reserves         = 0.01; %decent/approach/taxi



clmax = 1.8;

clmin = 0.2; 

cd0 = cdo;% with payload drag coefficient

cd0c = 0.025; % no external payload drag coefficient



%Mission Req

% requirement 1 MAX PL and 40 hrs En

mission_1_x = [40 40];

mission_1_y = [0 60e3];



% requirement 2 half PL and 50 hrs En

mission_2_x = [50 50];

mission_2_y = [0 60e3];


% requirement 3 rotor takeoff

mission_3_x = [40 55];

mission_3_y = [23e3 23e3];



%Vectors of PL and Endurance

pl_vector = linspace(500,3500,4);

en_vector = linspace(40,55,10); 

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

assumptions = sprintf('cdo: %.4f, TSFC: %.6f, e: %.2f, alt: %.0f ft, S: %.0f m^2, AR: %.0f, Cruise: %.0f kts , Loiter: %.0f kts ', cdo,TSFC,e,cruise_alt*3.281,S,AR,convvel(v_cruise,'m/s','kts'),convvel(v_loiter,'m/s','kts'));

annotation('textbox',...
    [0.33 0.86 0.21 0.05],...
    'String',assumptions,...
    'LineStyle','none',...
    'FontSize',12,...
    'FitBoxToText','off');



% Plot mission lines

plot(mission_1_x,mission_1_y,'--','color','k','linewidth',3)

plot(mission_2_x,mission_2_y,'--','color','k','linewidth',3)

plot(mission_3_x,mission_3_y,'--','color','k','linewidth',3)

tow = TOW_results(~isinf(TOW_results));

ylim([min(min(tow))-1e4,max(max(tow))+1e4])

xlim([40,55])

legend(legend_labels,'location','NW')
