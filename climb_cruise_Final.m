%Enviromental

    g    = 9.81;

    T0   = 288;      % +30; ISA + 30

    L_r  = 0.0065; %lapse rate for temp

    R    = 287;

    rho0 = 1.225;    %1.10949; ISA + 30

    P0   = 101300;


%Initalise variables

    time_res = 1/60; % 1 min res


%Set inital conditions (starting conditions)
    clear('w','t','d','alt','v','ROC_max')

    w(1)       = TOW;

    t(1)       = 0; % set time to 0  
    
    d(1)       = 0;

    alt(1)     = 0; % set altitude to 0
    
    clear fused;

    counter_2  = 0; % counter for sake of PL drop

    PL_dropped = 0;
    
    i = 1;
        
     while w(i) > (empty_weight + 700)

        %determine ambients (Linear - No account for Tropo)

        temp = T0-alt(i)*L_r;

        press= P0*(temp/T0)^(g/L_r/R);

        rho(i)  = press/(R*temp);
        
        % Determine if drop has occured
        if t(i) >= endurance_target/2  && counter_2 == 0 

            PL_dropped = PL - 500;

            w(i) = w(i) - PL_dropped; %loose PL

            counter_2 = counter_2 + 1;

        end
         
        
        % assumed operation is CL = 1.0 (an
        if mission_index < 3 
        
            cl(i) = min_drag_factor;
       
        else
            
            cl(i) = min_drag_factor;
            
        end
        

        L(i)    = w(i)*g; %N - Small angle assumption for climb

        if alt < cruise_alt
            
                v(i) = (L(i)/(0.5*rho(i)*S*cl(i)))^0.5;
                
                engine_characteristics_pmax
                
                % Calculate Cd
                cd_function
            
                D(i) = cd(i)*(0.5*rho(i)*S*v(i)^2);
                
                %determine max climb performance
                
                Thrust_max = P_max(i)*prop_n/v(i);
                ROC_max = v(i)*(Thrust_max - D(i))/L(i);

                if ROC_max > target_roc
                    
                    roc(i) = target_roc;
                    Thrust(i) = roc(i)*L(i)/v(i)+D(i);
                    Power(i) = Thrust(i)*v(i)./prop_n;
                 
                    
                else
                    
                    roc(i)    = ROC_max;
                    Thrust(i) = roc(i)*L(i)/v(i)+D(i);
                    Power(i)  = P_max(i);
                    
                end
                    
                alt(i+1) = alt(i) + (Thrust(i) - D(i))*v(i)/L(i)*time_res*3600;  %m, time res is in hrs
      
        else %Cruise or Loiter (aircraft needs to be at cruising level)
            
            
            v(i)= (L(i)/(0.5*rho(i)*S*cl(i)))^0.5;
            
         
            cd_function;
            
            Thrust(i) = cd(i)*(0.5*rho(i)*S*v(i)^2);
            
            cruising_prop_n

            Power(i) = Thrust(i)*v(i)/prop_n;
            
            alt(i+1) = alt(i);
            
        end
        
        engine_characteristics_BSFC
       
        w(i+1) = w(i) - Power(i)*BSFC(i)*time_res*3600; %end of segement fuel
        
        d(i+1) = d(i) + v(i)*time_res*3600;

        t(i+1) = t(i) + time_res; 
        
        i = i+1;

     end
    
     fused = sum(Power.*BSFC*time_res*3600);
     t_run = t;
     [ ~, loiter_start_t] = min(abs( d-loiter_point));
     loiter_start_t = t(loiter_start_t);
     
     [ ~, loiter_end_t ] = min(abs( d-(max(d)-loiter_point)));
     loiter_end_t = t(loiter_end_t);
      
     useful_loiter = loiter_end_t - loiter_start_t;
       