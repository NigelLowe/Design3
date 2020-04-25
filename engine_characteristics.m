a = (287*temp*1.4)^0.5;
mach(i)      = v(i)/a;

%find closest mach no
[ ~, mach_ii] = min(abs(FATE_mach-mach(i)));
[ ~, alt_ii]  = min(abs(FATE_alt-alt(i)));

%determine where next point is
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

v_P_max = [FATE_power(mach_ind(1),alt_ind(1)),FATE_power(mach_ind(1),alt_ind(2));...
            FATE_power(mach_ind(2),alt_ind(1)),FATE_power(mach_ind(2),alt_ind(2))];

v_BSFC = [FATE_BSFC(mach_ind(1),alt_ind(1)),FATE_BSFC(mach_ind(1),alt_ind(2));...
            FATE_BSFC(mach_ind(2),alt_ind(1)),FATE_BSFC(mach_ind(2),alt_ind(2))];
    
%Find in spreadsheet (iterpolate)

P_max(i) = interp2(alt_ref,mach_ref,v_P_max,alt(i),mach(i))*1e3;
BSFC(i)  = interp2(alt_ref,mach_ref,v_BSFC,alt(i),mach(i))/(60*60*1e3);


                
