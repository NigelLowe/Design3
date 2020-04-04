%% Call all Albatross Parameters - Maritime

%Wing + Aero Parameters

S   = 48.4;        %m2
AR  = 10;            
b   = sqrt(S*AR); %m
c   = S/b;

e   = 0.85;
cdo = 0.034;
k   = 1/(pi*e*AR);

%Engine
% TSFC = 14*1e-6;   %kg/(s.N)
prop_n = 0.8;

empty_weight = 6500;  % Used when BEW is fixed

%initalise deviations
isa_dev = 0;   %degrees kel
alt_dev = 0; %ft

%mission parameters
reach_toc       = 4;  %hrs
cruise_alt      = 35000/3.281; %m
loiter_point    = 1500e3; %4000 km
v_cruise        = convvel(250,'kts','m/s');
%v_loiter        = convvel(200,'kts','m/s');

target_roc      = cruise_alt/reach_toc/3600; %m/s
ld_climb        = 1/(2*sqrt(k*cdo)); %min l_d
cl_climb        = sqrt(pi*AR*cdo*e);

clmax = 1.8;
clmin = 0.2; 
cd0 = cdo;    % with payload drag coefficient
cd0c = 0.028; % no external payload drag coefficient

%Mission Req (for plotting)

% requirement 1 MAX PL

mission_1_x = [22 22];

mission_1_y = [0 18.0e3];

% requirement 2 50 hrs (reduced)

mission_2_x = [32 32];

mission_2_y = [0 18.0e3];

% requirement 3 rotor takeoff

mission_3_x = [10 55];

mission_3_y = [18.0e3 18.0e3];

x_lim_plot = [22-1 32+1];
