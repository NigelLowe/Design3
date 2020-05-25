%Calculate Cdo

if mission_index == 3 || mission_index == 4
   
    if counter_2 == 1
        
        cd(i) = 0.0345*cl(i).^2 - 0.0095*cl(i) + 0.0391;
        
    else
        
        cd(i) = 0.0294*cl(i).^2 + 0.0044*cl(i) + 0.04;
        
    end
    
elseif mission_index == 1 || mission_index == 2
    
    if counter_2 == 1
        
        cd(i) = 0.0294*cl(i).^2 - 0.0012*cl(i) + 0.055;
        
    else
        
        cd(i) = 0.022*cl(i).^2 + 0.0076*cl(i) + 0.057;
        
    end
    
end