%Enviromental
g    = 9.81;
T0   = 288;
L_r  = 0.0065; %lapse rate for temp
R    = 287;
rho0 = 1.225;
P0   = 101300;

%Initalise variables
TOW     = 30e3; %kg (guess starting weight)
TOW_new = 40e3;
n_time_pts = 3e4; 
time_res   = endurance_target/n_time_pts * 3600; %s

while abs(TOW_new - TOW)/TOW > 1e-5 
    w_initial = TOW;
    w = TOW;
    t_cruise = 0;
    t_total  = 0;
    alt      = 0;
    d_dist = 0;
        
    %for i = 1:n_time_pts+1
    while t_cruise < endurance_target*3600
        %determine ambients
        temp = T0-alt*L_r;
        press= P0*(temp/T0)^(g/L_r/R);
        rho  = press/(R*temp);
        
        L    = w*g; %N - Small angle assumption for climb
        
        % determine if aircraft in climb or cruise phase
        if alt < cruise_alt
            %Climb phase - Climb at min l_d
            D = L/ld_climb;
            v = (L/(0.5*cl_climb*rho*S))^0.5; %m/s
            
            %Determine Thrust req. (from excess power eqn.)
            Thrust = target_roc*L/v+D;
         
            %ROC
            roc = target_roc;
            %(roc = (Thrust - D)*v/L; %m/s)
            alt = alt + target_roc*time_res;  %m, time res is in hrs

        else %cruise at 80% Vmax
            v = v_cruise;
            roc = 0;
            cl   = L/(0.5*rho*S*v^2);
            cd   = cdo + k*(cl-clmin)^2;
            Thrust    = cd*(0.5*rho*S*v^2); % (T = D - dont need D variable later)
            if d_dist >= 4000e3    
                t_cruise = t_cruise + time_res; 
            end
        end

        %Fuel
        % (ff    = Thrust*TSFC)
        fused = Thrust*TSFC*time_res; %kg, as time res is in hrs

        w = w - fused; %end of segement fuel
        t_total = t_total + time_res;
        d_dist = d_dist + v*cos(roc)*time_res;
    end
    
    % weights    
    % Calc new TOW with mass fraction approach.
    % (empty_weight = 1.66*TOW^0.815) %As per texts
    % (w_initial - w) is fuel
    
    TOW_new = (1.66*TOW^0.815 + (w_initial - w) + PL)/(1-500/TOW);
    %TOW_new = (8700 + (w_initial - w) + PL)/(1-500/TOW); % Used when BEW is fixed
    
    % Take mean of 2 values for stability
    TOW = (TOW + TOW_new)/2;
    Power = rho0/rho * Thrust*v_cruise/0.85;
end