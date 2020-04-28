a = (287*temp*1.4)^0.5;
mach(i)      = v(i)/a;

%find closest mach no, altitude and power
[ ~, mach_ii] = min(abs(FATE_mach-mach(i)));
[ ~, alt_ii]  = min(abs(FATE_alt-alt(i)));

%determine where next point is (up or down?)
if (FATE_mach(mach_ii) - mach(i)) > 0
    
    mach_ref(1) = FATE_mach(mach_ii-1);
    mach_ref(2) = FATE_mach(mach_ii);
    
    mach_ind(1) = mach_ii-1;
    mach_ind(2) = mach_ii;
    
else
    
    mach_ref(1) = FATE_mach(mach_ii);
    mach_ref(2) = FATE_mach(mach_ii+1);   
    
    mach_ind(1) = mach_ii;
    mach_ind(2) = mach_ii+1;
    
end
  
%determine where next point is (left or right?)
if (FATE_alt(alt_ii) - alt(i)) > 0
    
    alt_ref(1) = FATE_alt(alt_ii-1);
    alt_ref(2) = FATE_alt(alt_ii);
    
    alt_ind(1) = alt_ii-1;
    alt_ind(2) = alt_ii;
 
else
    
    alt_ref(1) = FATE_alt(alt_ii);
    alt_ref(2) = FATE_alt(alt_ii+1);  
    
    alt_ind(1) = alt_ii;
    alt_ind(2) = alt_ii+1;
    
end    

%Return a vector of available throttle settings for given mach/alitude
%pair, transform from 3D to 2D
Power_vec = permute(FATE_Power(mach_ind(1),alt_ind(1),:),[3 2 1])*1e3;
[ ~, Power_ii]  = min(abs(Power_vec - Power(i)));

%Determine if nex
if (Power_vec(Power_ii) - Power(i)) > 0
    
    power_ref(1) = Power_vec(Power_ii-1);
    power_ref(2) = Power_vec(Power_ii);
    
    power_ind(1) = Power_ii-1;
    power_ind(2) = Power_ii;
 
else
    
    power_ref(1) = Power_vec(Power_ii);
    power_ref(2) = Power_vec(Power_ii+1);  
    
    power_ind(1) = Power_ii;
    power_ind(2) = Power_ii+1;
    
end  


v_BSFC(1:2,1:2,1) = [FATE_BSFC(mach_ind(1),alt_ind(1),power_ind(1)),FATE_BSFC(mach_ind(1),alt_ind(2),power_ind(1));...
                     FATE_BSFC(mach_ind(2),alt_ind(1),power_ind(1)),FATE_BSFC(mach_ind(2),alt_ind(2),power_ind(1))];

v_BSFC(1:2,1:2,2) = [FATE_BSFC(mach_ind(1),alt_ind(1),power_ind(2)),FATE_BSFC(mach_ind(1),alt_ind(2),power_ind(2));...
                     FATE_BSFC(mach_ind(2),alt_ind(1),power_ind(2)),FATE_BSFC(mach_ind(2),alt_ind(2),power_ind(2))];
                 
%Find in spreadsheet (iterpolate)

BSFC(i)  = interp3(alt_ref,mach_ref,power_ref,v_BSFC,alt(i),mach(i),Power(i))/(60*60*1e3);


                
