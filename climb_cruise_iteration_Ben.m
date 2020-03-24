%Enviromental
g    = 9.81;
T0   = 288;
L_r  = 0.0065; %lapse rate for temp
R    = 287;
rho0 = 1.225;
P0   = 101300;

%Initalise variables
TOW = 30e3; %kg (guess starting weight)
delta = 1;
n_time_pts = 1e5; 
delta_counter = 0; %to avoid infinte loop

while delta > 1e-5 && delta_counter < 1e4
    
    delta_counter = delta_counter + 1; 
    time_res   = endurance_target/n_time_pts; %hrs
    
    %Set inital conditions (starting conditions)
    w(1)       = TOW;
    i          = 0; % establish a counter for sake of referencing
    t          = 0; % set time to 0
    alt(1)     = 0; % set altitude to 0
    
    counter_2 = 0; %counter for sake of PL drop
    
    while t <= endurance_target
        
        i = i+1;
        
        %determine ambients (Linear - No account for Tropo)
        temp(i) = T0-alt(i)*L_r;
        press(i)= P0*(temp(i)/T0)^(g/L_r/R);
        rho(i)  = press(i)/(R*temp(i));
        
        time(i) = t; %hrs
        
        %determine if PL has been dropped? This only happens once
        if t >= endurance_target/2
            counter_2 = counter_2+1;
            
            if counter_2 == 1;
                w(i) = w(i) - 0.5*PL; %loose half PL
                cdo  = 0.025;  % drag improvement from external payload
            end
            
        end
        
        L(i)    = w(i)*g; %N - Small angle assumption for climb
        
        % determine if aircraft in climb or cruise phase
        if alt(i) < cruise_alt
            cl(i) = cl_climb;
            
            %Climb phase - Climb at min l_d
            D(i) = L(i)/ld_climb;
            v(i) = (L(i)/(0.5*cl(i)*rho(i)*S))^0.5; %m/s
            
            %Determine Thrust req. (from excess power eqn.)
            Thrust(i) = target_roc*w(i)*g/v(i)+D(i);
            
            %ROC
            roc(i) = (Thrust(i) - D(i))*v(i)/(w(i)*g); %m/s
            alt(i+1) = alt(i) + roc(i)*time_res*3600;  %m, time res is in hrs
            
            %Fuel
            ff(i)    = Thrust(i)*TSFC;
            fused(i) = ff(i)*time_res*3600; %kg, as time res is in hrs
            
            w(i+1) = w(i) - fused(i); %end of segement fuel
         
        else %cruise at 80% Vmax

            v(i) = 102.8889/0.8*0.9; %200 kts
            cl(i)= L(i)/(0.5*rho(i)*S*v(i)^2);
            
            cd   = cdo + k*cl(i)^2;
            D(i) = cd*(0.5*rho(i)*S*v(i)^2);
            
            L_D(i) = L(i)/D(i);
            
            Thrust(i) = D(i); %N
            
            %ROC
            roc(i) = (Thrust(i) - D(i))*v(i)/(w(i)*9.81);
            alt(i+1) = alt(i) + roc(i)*time_res*3600;
            
            %Fuel
            ff(i)    = Thrust(i)*TSFC;
            fused(i) = ff(i)*time_res*3600;   %kg, as time res is in hrs
          
            w(i+1) = w(i) - fused(i); %end of segement fuel
        
        end
    
        t = t + time_res; 
        
    end
    
    % weights
    %empty_weight = 1.66*TOW^0.815; %As per texts
    empty_weight = 6000;            % Used when BEW is fixed
    
    fuel         = sum(fused); %account for PL loss too;
    
    % Calc new TOW with mass fraction approach.
    TOW_new = (empty_weight + fuel + PL)/(1-500/TOW);
    
    % Calculate new delta for while critera
    delta = abs(TOW_new - TOW)/TOW;
    
    % Take mean of 2 values for stability
    TOW = (TOW + TOW_new)/2;
    
end
