clear all, 
close all,
clc

% DESIGN 3 - Performance Script - Group 3

% Call Parameters
albatross_parameters_maritime
%albatross_parameters_airfield

%Vectors of PL and Endurance

pl_vector = linspace(500,500,1);

en_vector = linspace(32,32,1); 

%Wait bar for sanity

hwait = waitbar(0,'Please wait...');

h_bar = 1;

for m_index = 1:length(pl_vector)

    

    PL = pl_vector(m_index);

    for n_index = 1:length(en_vector)

       
        endurance_target = en_vector(n_index);

        climb_cruise_iteration

        TOW_results(m_index,n_index) = TOW;
        BEW_results(m_index,n_index) = empty_weight;
        clmax_results(m_index,n_index)  = max(cl);
        T_max_results(m_index,n_index)  = max(Thrust);
        P_max_results(m_index,n_index)  = max((Thrust.*v)/prop_n);
        WL_max_results(m_index,n_index) = TOW/S;
      
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

xlim(x_lim_plot)

legend(legend_labels,'location','NW')
