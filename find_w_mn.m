% aerodynamics assignment 2 
% question 4 
% 
% this function finds the w_mn matrix for the VLM of a wing 
% 
% inputs: 
%           N - number of panels
%           x_1n - vector of x coordiantes of n1 points 
%           y_1n - vector of y coordinates of n1 points 
%           x_2n - vector of x coordiantes of n2 points  
%           y_2n - vector of y coordiantes of n2 points 
%           x_m - vector of x coordinates of m points
%           y_m - vector of y coordinates of m points
% 
% outputs: 
%           w_mn matrix 


function w_mn = find_w_mn(N, x_1n, y_1n, x_2n, y_2n, x_m, y_m)

    % initialise the w_{mn} matrix 
    w_mn = zeros(N,N); 

    % loop through each panel (rows first, then columns) 
    for index_N = 1:N

        % loop through each panel again to find influence of n on m
        for index_N_nested = 1:N   

            A = 1/((x_m(index_N) - x_1n(index_N_nested)) * (y_m(index_N) - y_2n(index_N_nested)) - (x_m(index_N) - x_2n(index_N_nested)) * (y_m(index_N) - y_1n(index_N_nested))); 

            B = ((x_2n(index_N_nested) - x_1n(index_N_nested)) * (x_m(index_N) - x_1n(index_N_nested)) + (y_2n(index_N_nested) - y_1n(index_N_nested)) * (y_m(index_N) - y_1n(index_N_nested)))/...
                sqrt((x_m(index_N) - x_1n(index_N_nested))^2 + (y_m(index_N) - y_1n(index_N_nested))^2); 

            C = - ((x_2n(index_N_nested) - x_1n(index_N_nested)) * (x_m(index_N) - x_2n(index_N_nested)) + (y_2n(index_N_nested) - y_1n(index_N_nested)) * (y_m(index_N) - y_2n(index_N_nested)))/...
                sqrt((x_m(index_N) - x_2n(index_N_nested))^2 + (y_m(index_N) - y_2n(index_N_nested))^2);

            D = 1/(y_1n(index_N_nested) - y_m(index_N)) * (1 + (x_m(index_N) - x_1n(index_N_nested))/sqrt((x_m(index_N) - x_1n(index_N_nested))^2 + (y_m(index_N) - y_1n(index_N_nested))^2)); 

            E = - 1/(y_2n(index_N_nested) - y_m(index_N)) * (1 + (x_m(index_N) - x_2n(index_N_nested))/sqrt((x_m(index_N) - x_2n(index_N_nested))^2 + (y_m(index_N) - y_2n(index_N_nested))^2)); 

            w_mn(index_N,index_N_nested) = (1/(4*pi)) * (A * (B+C) + D + E);

        end 
    end

end 

