function [U_Lon, U_Lat] = Controls(t, mnvre)

% Longitudinal control
U_Lon = zeros(2,1);
% Lateral-directional control
U_Lat = zeros(2,1);
% Manoeuvre duration
    if (t<1)||(t>1.5)
        % Longitudinal control
        U_Lon = zeros(2,1);
        % Lateral-directional control
        U_Lat = zeros(2,1);
            % Elevator deflection 
            elseif mnvre == 1
                U_Lon (2) = deg2rad(5) ;
            % Aileron deflection
            elseif mnvre == 2
                U_Lat (1) = deg2rad(5) ;
            % Rudder deflection
            elseif mnvre == 3
                U_Lat (2) = deg2rad(5) ;
    end
end

