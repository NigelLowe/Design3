if ~exist('plotOtherGraphs','var') % if statement for this file to be used in other functions
    clear all
end 
close all,
clc

% DESIGN 3 - Performance Script - Group 3
load('FATE_parameters.mat')

%Define Mission Constants
mission_PL           = [3500,500,3500,500];
mission_weight       = [18.5e3,18.5e3,22.5e3,22.5e3];
endurance_target_mat = [25,36,50,60];
cruising_alt_target  = 40000;

missionIndecies = 1:4;
if exist('plotOtherGraphs','var') 
    mission_PL = pl_num;
    mission_weight = tow_num;
    endurance_target_mat = en_num;
    missionIndecies = 1;
end

for mission_index = missionIndecies
   
        PL               = mission_PL(mission_index);
        TOW              = mission_weight(mission_index);
        endurance_target = endurance_target_mat(mission_index);
        
        %Call relevant Parameters
        if ~exist('plotOtherGraphs','var')
            if mission_index < 3
                albatross_parameters_maritime
            else
                albatross_parameters_airfield
            end
        end
        
        %Refine Cruising Alt Sweep
        cruising_alt_mat = convlength(linspace(30e3,40e3,10),'ft','m');
        
        % Min Drag Speed factor Sweep
        min_drag_factor_mat = linspace(0.5,1.5,10);
        
        for m_index = 1:length(cruising_alt_mat)
            
            for n_index = 1:length(min_drag_factor_mat)
                
                cruise_alt      = cruising_alt_mat(m_index);
                min_drag_factor = min_drag_factor_mat(n_index);
             
        climb_cruise_CDR

        results_TOW(m_index,n_index,mission_index)               = TOW;
        results_endurance(m_index,n_index,mission_index)         = max(t);
        results_range(m_index,n_index,mission_index)             = max(d);
        results_useful_loiter(m_index,n_index,mission_index)     = useful_loiter;
        results_average_velocity(m_index,n_index,mission_index)  = mean(v);
        results_max_cl(m_index,n_index,mission_index)            = max(cl);
        results_mean_weight(m_index,n_index,mission_index)       = mean(w);
        
            end
        
        end
        
 %create ALT vs Vel Contour (with endurance)
    if ~exist('plotOtherGraphs','var')
         figure(mission_index)
         hold on

            [~,alt_grid] = meshgrid(cruising_alt_mat,cruising_alt_mat);
            [c1,h1] = contour(results_average_velocity(:,:,mission_index)*1.9438,alt_grid*3.2808,results_endurance(:,:,mission_index),'LineWidth',2);
            clabel(c1,h1)

            %calculate/plot stall speeds
            for stall_index = 1:length(cruising_alt_mat)

                [~, ~, ~, rho] = atmosisa(cruising_alt_mat(stall_index));
                stall_speed(stall_index) = ((mean(results_mean_weight(stall_index,:,mission_index))*9.81)/(0.5*rho*clmax*S))^0.5;

            end

           plot(stall_speed*1.9438,cruising_alt_mat*3.2808,'--r','LineWidth',2)

            %calculate/plot mach 0.5

            for mach_index = 1:length(cruising_alt_mat)

                [~, a, ~, ~] = atmosisa(cruising_alt_mat(mach_index));
                 mach_limit_speed(mach_index) = 0.5*a;

            end

            plot(mach_limit_speed*1.9438,cruising_alt_mat*3.2808,'--b','LineWidth',2)

            %locate design point (alt and velocity)
            [ ~, alt_ii] = min(abs(cruising_alt_mat-cruising_alt_target/3.2808));
            [~, vel_ii] = min(abs(results_endurance(alt_ii,:,mission_index)-endurance_target_mat(mission_index)));

            %Plot design point
            scatter([results_average_velocity(alt_ii,vel_ii,mission_index)*1.9438],[cruising_alt_target],'MarkerFaceColor',[0 0 0],'Marker','pentagram','SizeData',1000)

            %format figure
            set(gca,'fontsize',14)
            xlabel('True Airspeed (kts.)')
            ylabel('Altitude (ft.)')
            grid on

            ylim([(min(cruising_alt_mat*3.2808) - 1e3) (max(cruising_alt_mat*3.2808) + 1e3)]);

            %add legend
            legend('Endurance (hrs)','mean stall line','Mach 0.5','Design Point','location','best')


     %create ALT vs Vel Contour (with range)

         figure(mission_index+4)
         hold on
            [c2,h2] = contour(results_average_velocity(:,:,mission_index)*1.9438,alt_grid*3.2808,results_range(:,:,mission_index)*1e-3,'LineWidth',2);
            clabel(c2,h2)

        % Plot stall line
        plot(stall_speed*1.9438,cruising_alt_mat*3.2808,'--r','LineWidth',2)


        % Plot high speed line
        plot(mach_limit_speed*1.9438,cruising_alt_mat*3.2808,'--b','LineWidth',2)

            scatter([results_average_velocity(alt_ii,vel_ii,mission_index)*1.9438],[cruising_alt_target],'MarkerFaceColor',[0 0 0],'Marker','pentagram','SizeData',1000)

            set(gca,'fontsize',14)
            xlabel('True Airspeed (kts.)')
            ylabel('Altitude (ft.)')
            grid on

            ylim([(min(cruising_alt_mat*3.2808) - 1e3) (max(cruising_alt_mat*3.2808) + 1e3)]);
            legend('Air Range (km.)','mean stall line','Mach 0.5','Design Point','location','best')
    end
  
    if mission_index == 1
        saveas(figure(1),'Mission_1a_Endurance','png')
        saveas(figure(5),'Mission_1a_Range','png')
        
    elseif mission_index == 2
        
        saveas(figure(2),'Mission_1b_Endurance','png')
        saveas(figure(6),'Mission_1b_Range','png')
        
    elseif mission_index == 3
        
        saveas(figure(3),'Mission_2a_Endurance','png')
        saveas(figure(7),'Mission_2a_Range','png')
        
    elseif mission_index == 4
        saveas(figure(4),'Mission_2b_Endurance','png')
        saveas(figure(8),'Mission_2b_Range','png')
       
    end
    
end
