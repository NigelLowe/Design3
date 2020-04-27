% aerodynamics assignment 2 
% question 4 
% vortex lattice method for a rectangular wing 

clear;


set(0, 'defaultAxesFontSize', 20)
set(0, 'defaultlinelinewidth', 3);
set(0,'DefaultAxesXGrid','on', 'DefaultAxesYGrid','on', 'DefaultAxesZGrid','on')


%% GIVE THIS DATA TO CREATE THE WING YOU WANT 

% external dimensions of rectangular wing 
c_bar = 6;        % wing chord in metres 
b = 30;             % wing span in metres 
s = b/2;            % half wing span in metres 

% set number of rows to discretise the wing into 
total_rows = 20; 

% set number of columns to discretise the wing into 
total_cols = 25;

% set the taper ratio (c_t/c_r)
lambda = 0.75; 

% set the sweep angle
sweep_angle = deg2rad(10); 

% set the angle of attack of the wing 
% NOTE: this is restricted to 0-2 degrees
alpha = deg2rad(5); 

% set the freestream velocity in m/s 
U = 50; 

% density of air in kg/m^3
rho = 1.225; 

% set the dihedral angle 
gamma = deg2rad(4);


%% DONT NEED TO CHANGE THIS STUFF 

% number of panels 
N = total_rows * total_cols; 

% find the chordwise length of one panel 
c_bar_N = c_bar/total_rows; 

% find the spanwise length of one panel 
s_N = s/total_cols; 

% Initialise the vectors that contain all the 1n points on the wing (x and y coordiantes). 
% x and y coordiantes separately 
x_1n = zeros(N,1);
y_1n = zeros(N,1); 

% Initialise the vectors that contain all the 2n points 
% x and y coordinates separately 
y_2n = zeros(N,1); 

% Initialise the vectors that contain all the m points 
% x and y coordinates separately 
x_m = zeros(N,1); 
y_m = zeros(N,1); 

% Note: the origin is at the root leading edge, with x axis along the chord
% and y axis along the span. 

% Note: panel numbering ---> rows first and then columns 

%% populating the n1, n2 and m matrices 

% panel index number 
index_N = 0; 

% initialise the increment for the x coordinate for n1 
x_increment_n1 = 0.25 * c_bar_N;

% initialise the increment for the x coordinate for m 
x_increment_m = 0.75 * c_bar_N; 

% loop through each panel (rows first, then columns) 
for r = 1:total_rows
    
    % initialise the increment for the y coordinate for n1 
    y_increment_n1 = 0; 
    
    % initialise the increment for the y coordinate for n2
    y_increment_n2 = s_N; 
    
    % initialise the increment for the y coordinate for m 
    y_increment_m = 0.5 * s_N; 
    
    for c = 1:total_cols 

        % now we are at one panel so update the panel index number 
        index_N = index_N + 1;
        
        % x coordinate for the n1 vector 
        x_1n(index_N) = x_increment_n1; 
        
        % y coordiante for the n1 vector
        y_1n(index_N) = y_increment_n1; 
        
        % y coordinate for the n2 vector
        y_2n(index_N) = y_increment_n2; 
        
        % x coordinate for the m vector 
        x_m(index_N) = x_increment_m; 
        
        % y coordinate for the m vector 
        y_m(index_N) = y_increment_m; 
        
        % update the increment for the n1 y coordiante because we are going
        % into the next column 
        y_increment_n1 = y_increment_n1 + s_N; 
        
        % update the increment for the n2 y coord 
        y_increment_n2 = y_increment_n2 + s_N; 
        
        % update the increment for the m y coord 
        y_increment_m = y_increment_m + s_N; 
        
        
    end 
    
    % update the increment for the x coordinate because we are going into
    % the next row 
    x_increment_n1 = x_increment_n1 + c_bar_N; 
    
    % update the increment for the m x coord 
    x_increment_m = x_increment_m + c_bar_N; 
    
end 
        

% the x coordinates of the n2 points are equal to the x coordinates of the
% n1 points 
x_2n = x_1n;

% find the w_mn matrix 
w_mn = find_w_mn(N, x_1n, y_1n, x_2n, y_2n, x_m, y_m);

% this vector stores the spanwise length of each panels
% note that this does not change with taper and sweep, because the y
% coordiantes arent affected 
dy = ones(N,1) * s_N; 

% aoa of each panel 
alpha_bar = ones(N,1) * alpha;

% find the circulation and the lift 
[cap_gamma, dL_rec_wing] = find_lift(w_mn, U, alpha_bar, rho, dy); 

% find the total lift on the rectangular wing 
L_rec_wing  = sum(dL_rec_wing);

% find the lift coefficient 
CL_rec_wing = L_rec_wing/(0.5 * rho * U^2 * c_bar * s); 

%% plot the points 

figure(1); 
clf; 
clf reset; 
plot(x_1n, y_1n, 'ro', 'Linewidth', 1); 
hold on; 
plot(x_2n, y_2n, 'bo', 'Linewidth', 1); 
hold on; 
plot(x_m, y_m, 'ko', 'Linewidth', 1);
grid on;
xlim([0,c_bar]);
xlabel('distance along chord (m)'); 
ylabel('distance along span (m)'); 
title('Rec wing (origin at root LE)');
legend('n1', 'n2', 'm'); 

%% plot the lift distribution on the rec wing 

figure(4); 
clf; 
clf reset; 
plot3(x_1n, y_1n, dL_rec_wing, 'o'); 
xlabel('distance along chord (m)'); 
ylabel('distance along span (m)'); 
zlabel('lift (N)');
title('Lift Distribution for rectangular wing'); 
grid on; 

%% plotting panels properly 

x_rec_plot = (0:c_bar_N:c_bar)'; 
y_rec_plot = (0:s_N:s)';




%% adding taper 

% note: taper ratio goes from 1 to lambda as you move from root to tip

% all spanwise coordinates from root to tip
span_coord = 0: 0.5 * s_N :s; 

% all lambda values from root to tip 
lambda_vec = interp1([0, s], [1, lambda], span_coord); 

% all spanwise coordinates for n1, n2 and m separately
span_coord_n1 = 0:s_N:s-s_N;
span_coord_n2 = s_N:s_N:s;
span_coord_m = (s_N/2):s_N:(s-(s_N/2));

% find the lambda values for n1, n2 and m
lambda_vec_n1_row1 = interp1([0, s-s_N], [1, lambda_vec(end-2)], span_coord_n1)';
lambda_vec_n2_row1 = interp1([s_N, s], [lambda_vec(3), lambda], span_coord_n2)';
lambda_vec_m_row1 = interp1([s_N/2, (s-(s_N/2))], [lambda_vec(2), lambda_vec(end-1)], span_coord_m)';

% initialise the lambda vectors for n1, n2 and m for all the rows 
lambda_vec_n1 = []; 
lambda_vec_n2 = []; 
lambda_vec_m = []; 

% loop through the number of rows minus 1
% concatenation 
for r = 1:total_rows
    lambda_vec_n1 = [lambda_vec_n1;lambda_vec_n1_row1];
    lambda_vec_n2 = [lambda_vec_n2;lambda_vec_n2_row1];
    lambda_vec_m = [lambda_vec_m;lambda_vec_m_row1];
    
end 

% y coordinate stays the same with taper 
% only the x coordinate changes i.e. chordwise positions 

% update the x coordinates to account for taper 
x_1n_tap = x_1n .* lambda_vec_n1;
x_2n_tap = x_2n .* lambda_vec_n2;
x_m_tap = x_m .* lambda_vec_m;

%% plot the results for taper 

figure(2); 
clf; 
clf reset; 
plot(x_1n_tap, y_1n, 'ro', 'Linewidth', 1); 
hold on; 
plot(x_2n_tap, y_2n, 'bo', 'Linewidth', 1); 
hold on; 
plot(x_m_tap, y_m, 'ko', 'Linewidth', 1);
grid on;
xlim([0,c_bar]);
xlabel('distance along chord (m)'); 
ylabel('distance along span (m)'); 
title('Tapered wing (origin at root LE)');
legend('n1', 'n2', 'm'); 

%% adding sweep 

% again y coordinate stays the same with taper 
% only the x coordinate changes i.e. chordwise positions 

% update the x coordinates to account for sweep angle
x_1n_sweep = x_1n_tap + y_1n * tan(sweep_angle); 
x_2n_sweep = x_2n_tap + y_2n * tan(sweep_angle);
x_m_sweep = x_m_tap + y_m * tan(sweep_angle);


%% plot the results for sweep

figure(3); 
clf; 
clf reset; 
plot(x_1n_sweep, y_1n, 'ro', 'Linewidth', 1); 
hold on; 
plot(x_2n_sweep, y_2n, 'bo', 'Linewidth', 1); 
hold on; 
plot(x_m_sweep, y_m, 'ko', 'Linewidth', 1);
grid on;
xlim([0,c_bar]);
xlabel('distance along chord (m)'); 
ylabel('distance along span (m)'); 
title('Tapered and Swept back wing (origin at root LE)');
legend('n1', 'n2', 'm'); 



%% find the lift for the tapered, swept back wing 

% aoa of each panel 
alpha_bar = ones(N,1) * alpha;

% find the w_mn matrix 
w_mn = find_w_mn(N, x_1n_sweep, y_1n, x_2n_sweep, y_2n, x_m_sweep, y_m);

% find the circulation and the lift 
[cap_gamma, dL_ts_wing] = find_lift(w_mn, U, alpha_bar, rho, dy); 

% find the total lift on the tapered, swept back wing
L_ts_wing = sum(dL_ts_wing); 


%% plot the lift distribution on the tapered, swept back wing 

figure(5); 
clf; 
clf reset; 
plot3(x_1n_sweep, y_1n, dL_ts_wing, 'o'); 
xlabel('distance along chord (m)'); 
ylabel('distance along span (m)'); 
zlabel('lift (N)'); 
title('Lift Distribution for wing with taper, sweep'); 
grid on; 

%% trying to plot a surface plot of x,y,z with lift as a colormap

% dL_rec_wing is a Nx1 vector, which gives the lift on each panel 

% Create new arrays for surface plotting
L = zeros(2);
z_surf = zeros(size(x_1n));

% Initialise index for lift vector element

VLM_surface_plotter_3D(c_bar, ones(1,2*total_cols+1), x_1n, x_2n, y_rec_plot, z_surf, dL_rec_wing, total_rows, total_cols)
title('Lift Distribution for rectangular wing'); 


%% Plot lift as a surface plot in a new figure
% Initialise index for lift vector element

% Get an array for the chord at each y node
VLM_surface_plotter_3D(c_bar, lambda_vec, x_1n_sweep, x_2n_sweep, y_rec_plot, z_surf, dL_ts_wing, total_rows, total_cols)
title('Lift Distribution for tapered swept wing'); 
