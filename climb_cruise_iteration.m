%Enviromental
g    = 9.81;
T0   = 288;
L_r  = 0.0065; %lapse rate for temp
R    = 287;
rho0 = 1.225;
P0   = 101300;

%Initalise variables
TOW     = endurance_target/2 * 1e3; %kg (guess starting weight)
TOW_new = 2*TOW;
% delta = 1;
n_time_pts = 1e4; 
delta_counter = 0; %to avoid infinte loop
time_res   = endurance_target/n_time_pts; %hrs

while abs(TOW_new - TOW)/TOW > 1e-5
    
    delta_counter = delta_counter + 1; 
    
    %Set inital conditions (starting conditions)
    w          = TOW;
    w_initial  = TOW;
    t          = 0; % set time to 0
    alt        = 0; % set altitude to 0
    counter_2  = 0; %counter for sake of PL drop
    PL_dropped = 0;
    
    while t <= endurance_target
  
        %determine ambients (Linear - No account for Tropo)
        temp = T0-alt*L_r;
        press= P0*(temp/T0)^(g/L_r/R);
        rho  = press/(R*temp);

        if t >= endurance_target/2  && counter_2 == 0 
            PL_dropped = 0.5*PL;
            w = w - PL_dropped; %loose PL
            cdo  = cd0c;  % drag improvement from external payload
            counter_2 = counter_2 + 1;
        end 
        
        L    = w*g; %N - Small angle assumption for climb
        
        % determine if aircraft in climb or cruise phase
        if alt < cruise_alt
            cl = cl_climb;
            
            %Climb phase - Climb at min l_d
            D = L/ld_climb;
            v = (L/(0.5*cl*rho*S))^0.5; %m/s
            
            %Determine Thrust req. (from excess power eqn.)
            Thrust = target_roc*w*g/v+D;
            
            alt = alt + (Thrust - D)*v/L*time_res*3600;  %m, time res is in hrs
            
        else %cruise at 80% Vmax

            cl= L/(0.5*rho*S*v_cruise^2);
            cd   = cdo + k*(cl-clmin)^2;
            Thrust = cd*(0.5*rho*S*v_cruise^2);
        
        end
    
        %Fuel
        %(ff = Thrust*TSFC; % fuel flow)
        %(fused(i) = Thrust*TSFC*time_res*3600; %kg, as time res is in hrs)

        w = w - Thrust*TSFC*time_res*3600; %end of segement fuel
            
        t = t + time_res; 
    end
    
    % weights
    empty_weight = 1.66*TOW^0.815; %As per texts
    %empty_weight = 6000;            % Used when BEW is fixed
    
    % Calc new TOW with mass fraction approach.
    TOW_new = (empty_weight + (w_initial - w) + (PL - PL_dropped))/(1-500/TOW);
    
    
    % Take mean of 2 values for stability
    TOW = (TOW + TOW_new)/2;
    
    cdo = cd0; % reset cdo foor next loop
    if delta_counter > 500 % stop computing if TOW gets too big
        TOW = Inf;
        break;
    end
end
