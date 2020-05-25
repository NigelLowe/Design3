% CALCUATE cruising PROP n through reference of airspeed and thrust req

%find closest thrust
[ ~, Thrust_ii] = min(abs(PROP_Thrust-Thrust(i)));
[ ~, v_ii]  = min(abs(PROP_airspeed-v(i)));

if (PROP_Thrust(Thrust_ii) - Thrust(i)) > 0
    
    Thrust_ref(1) = PROP_Thrust(Thrust_ii-1);
    Thrust_ref(2) = PROP_Thrust(Thrust_ii);
   
    Thrust_ind(1) = Thrust_ii-1;
    Thrust_ind(2) = Thrust_ii;
    
else
    
    Thrust_ref(1) = PROP_Thrust(Thrust_ii);
    Thrust_ref(2) = PROP_Thrust(Thrust_ii+1);  
    
    Thrust_ind(1) = Thrust_ii;
    Thrust_ind(2) = Thrust_ii+1;
    
end 

if (PROP_airspeed(v_ii) - v(i)) > 0
    
    v_ref(1) = PROP_airspeed(v_ii-1);
    v_ref(2) = PROP_airspeed(v_ii);
   
    v_ind(1) = v_ii-1;
    v_ind(2) = v_ii;
    
else
    
    v_ref(1) = PROP_airspeed(v_ii);
    v_ref(2) = PROP_airspeed(v_ii+1);  
    
    v_ind(1) = v_ii;
    v_ind(2) = v_ii+1;
    
end 

v_prop_n = [PROP_n_range(v_ind(1),Thrust_ind(1)),PROP_n_range(v_ind(1),Thrust_ind(2));...
           PROP_n_range(v_ind(2),Thrust_ind(1)),PROP_n_range(v_ind(2),Thrust_ind(2))];
       
prop_n = interp2(Thrust_ref,v_ref,v_prop_n,Thrust(i),v(i));

           
