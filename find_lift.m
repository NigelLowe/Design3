% aerodynamics assignment 2 
% question 4 
% 
% this function returns the circulation and lift for the wing 
% 
% inputs: 
%           w_mn        - influence matrix 
%           U           - freestream velocity in m/s
%           alpha_bar   - panel angles of attack (include camber effects) 
%           rho         - density in kg/m^3
%           dy          - spanwise length of each panel (Nx1)
%
% 
% outputs: 
%           cap_gamma	- vector of circulation for each panel (Nx1)
%           dL          - vector of lift acting on each panel (Nx1)
% 
% 

function [cap_gamma, dL] = find_lift(w_mn, U, alpha_bar, rho, dy) 
    
    % find the circulation vector
    cap_gamma = - inv(w_mn) * U * sin(alpha_bar); 
    
    % find the lift on each panel using Kutta Joukowski theorem
    % note we need to multiply by dy to get depth (3D effect)
    dL = cap_gamma * rho * U .* dy; 

end 








