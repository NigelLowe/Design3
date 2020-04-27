% Aerodynamics 1
% Assignment 2 Question 4 (VLM)
% Author: Jack Knight
%
% Function that generates a surface plot of the lift generated in each
% panel for Vortex Lattice Method (3D)
function VLM_surface_plotter_3D(c_bar, lambda_vec, x_1n, x_2n, y, z_1n, dL, rows, cols)

    
    % Get ready to plot on a new figure
    figure; 
    clf; 
    clf reset; 
    hold on;
    
    % Initialise an index for the current panel we are looking at
    index = 1;
    
    % Get the local chord at each y node
    lambda_nodes = lambda_vec(1:2:end);
    c_nodes = c_bar * lambda_nodes;

    % Get separate plots for each panel
    for i = 1:rows %x
        for j = 1:cols %y

            % Get the current chord length
            chord1 = c_nodes(j) / rows;
            chord2 = c_nodes(j+1) / rows;

            % Get matrixes for the x y and L for the current panel
            x_surf = [x_1n(index) - chord1/4, x_2n(index) - chord2/4;...
                      x_1n(index) + 3*chord1/4, x_2n(index) + 3*chord2/4];

            y_surf = [y(j:j+1)'; y(j:j+1)'];
            
            z_surf = [z_1n(index), z_1n(index) + (z_1n(2)-z_1n(1));...
                      z_1n(index), z_1n(index) + (z_1n(2)-z_1n(1))];
                  
            L = dL(index) * ones(2);

            % Plot the current panel
            surf(x_surf,y_surf,z_surf, L,'edgecolor','k','facecolor','interp');

            % Increment the count
            index = index + 1;
        end
    end
    
    % Add features to the plot
    view(0,90)
    col = colorbar;
    xlabel('Chord-wise Position x (m)')
    ylabel('Span-wise Position y (m)')
    zlabel('Height above Sea Level z (m)')
end