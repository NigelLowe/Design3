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

n_time_pts = 1e4; 

delta_counter = 0; %to avoid infinte loop

time_res   = endurance_target/n_time_pts; %hrs

% % Mission schedule (in hrs)
% t_loiter_start = (loiter_point/v_cruise)/3600;
% t_loiter_end   = (endurance_target - t_loiter_start);

while abs(TOW_new - TOW)/TOW > 1e-5

    delta_counter = delta_counter + 1; 

    %Set inital conditions (starting conditions)

    w          = TOW;

    w_initial  = TOW;

    t          = 0; % set time to 0
    
    d          = 0; %mark distance flown

    alt        = 0; % set altitude to 0

    counter_2  = 0; %counter for sake of PL drop

    PL_dropped = 0;
 
            
    while t <= endurance_target

        %determine ambients (Linear - No account for Tropo)

        temp = T0-alt*L_r;

        press= P0*(temp/T0)^(g/L_r/R);

        rho  = press/(R*temp);

        if t >= endurance_target/2  && counter_2 == 0 

            PL_dropped = PL - 500;

            w = w - PL_dropped; %loose PL

            cdo  = cd0c;  % drag improvement from external payload

            counter_2 = counter_2 + 1;

        end 
        
        %fly at min drag speed
        cl = sqrt(cdo/k);

        L    = w*g; %N - Small angle assumption for climb

        %P_max = (2E-05*(alt)^2 - 0.7208*(alt) + 8146.2)*1e3;

        if alt < cruise_alt

                v = (L/(0.5*rho*S*cl))^0.5;
                
                cd   = cdo + k*(cl-clmin)^2;
            
                D = cd*(0.5*rho*S*v^2);
                
                Thrust = target_roc*L/v+D;
                
                Power = Thrust*v/prop_n;
                
                alt = alt + (Thrust - D)*v/L*time_res*3600;  %m, time res is in hrs
      
        else %Cruise or Loiter (aircraft needs to be at cruising level)
            
            v= (L/(0.5*rho*S*cl))^0.5;
            
            cd   = cdo + k*(cl-clmin)^2;
            
            Thrust = cd*(0.5*rho*S*v^2);

            Power = Thrust*v/prop_n;
            
            alt = alt;
            
        end
        
        
        %Fuel

        %(ff = Thrust*TSFC; % fuel flow)

        %(fused(i) = Thrust*TSFC*time_res*3600; %kg, as time res is in hrs)
        
%         TSFC = 0.19/(3600*1e3)*v/prop_n*2;

        BSFC = 1.0*(1e-17*(alt)^4 - 1e-13*(alt)^3 + 1e-09*(alt)^2 + 2e-07*(alt) + 0.2285)*1/(60*60*1e3);

        w = w - Power*BSFC*time_res*3600; %end of segement fuel
        
        d = d + v*t;

        t = t + time_res; 

    end

    

    % weights

    %empty_weight = 1.2*1.66*TOW^0.815; %As per texts


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

    %Set inital conditions (starting conditions)

    w(1)       = TOW;

    t(1)       = 0; % set time to 0  
    
    d(1)       = 0;

    alt(1)     = 0; % set altitude to 0

    counter_2  = 0; %counter for sake of PL drop

    PL_dropped = 0;
    
    i = 1;

     while t(i) <= endurance_target

        %determine ambients (Linear - No account for Tropo)

        temp = T0-alt(i)*L_r;

        press= P0*(temp/T0)^(g/L_r/R);

        rho(i)  = press/(R*temp);

        if t(i)>=endurance_target/2  && counter_2 == 0 

            PL_dropped = PL - 500;

            w(i) = w(i) - PL_dropped; %loose PL

            cdo  = cd0c;  % drag improvement from external payload

            counter_2 = counter_2 + 1;

        end 
        
        %fly at min drag speed
        cl(i) = sqrt(cdo/k);

        L(i)    = w(i)*g; %N - Small angle assumption for climb

%         P_max = (2E-05*(alt(i))^2 - 0.7208*(alt(i)) + 8146.2)*1e3;

       if alt(i) < cruise_alt

                v(i) = (L(i)/(0.5*rho(i)*S*cl(i)))^0.5;
                
                cd(i)   = cdo + k*(cl(i)-clmin)^2;
            
                D(i) = cd(i)*(0.5*rho(i)*S*v(i)^2);
                
                Thrust(i) = target_roc*L(i)/v(i)+D(i);
                
                Power(i) = Thrust(i)*v(i)/prop_n;
                
                alt(i+1) = alt(i) + (Thrust(i) - D(i))*v(i)/L(i)*time_res*3600;  %m, time res is in hrs
      
        else %Cruise or Loiter (aircraft needs to be at cruising level)
            
            v(i)= (L(i)/(0.5*rho(i)*S*cl(i)))^0.5;
           
            cd(i)   = cdo + k*(cl(i)-clmin)^2;
            
            Thrust(i) = cd(i)*(0.5*rho(i)*S*v(i)^2);

            Power(i) = Thrust(i)*v(i)/prop_n;
            
            alt(i+1) = alt(i);
            
        end
        
        
        %Fuel

        %(ff = Thrust*TSFC; % fuel flow)

        
%         TSFC = 0.19/(3600*1e3)*v/prop_n*2;

        BSFC(i) = 1.0*(1e-17*(alt(i))^4 - 1e-13*(alt(i))^3 + 1e-09*(alt(i))^2 + 2e-07*(alt(i)) + 0.2285)*1/(60*60*1e3);
        
        fused(i) = BSFC(i)*Power(i)*time_res*3600; %kg, as time res is in hrs)
        
        TSFC(i) = BSFC(i)*Power(i)/Thrust(i);

        fused(i) = Thrust(i)*TSFC(i)*time_res*3600; %kg, as time res is in hrs)
        
        w(i+1) = w(i) - BSFC(i)*Power(i)*time_res*3600; %end of segement fuel
       
        t(i+1) = t(i) + time_res; 
        
        d(i+1) = d(i) + v(i)*time_res*3600;
        
        %advance referencing parameter
        i = i + 1; 

     end
     
     [ ~, loiter_start_t ] = min(abs( d-loiter_point));
     loiter_start_t = t(loiter_start_t);
     
     [ ~, loiter_end_t ] = min(abs( d-(max(d)-loiter_point)));
     loiter_end_t = t(loiter_end_t);
      
     useful_loiter = loiter_end_t - loiter_start_t;
       
