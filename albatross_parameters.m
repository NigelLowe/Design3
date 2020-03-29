%% Call all Albatross Parameters

%Wing Geometry

S   = 68;        %m2
AR  = 18;            
b   = sqrt(S*AR); %m
c   = S/b;

% Aerodynamic Parameters
e   = 0.85;
k   = 1/(pi*e*AR);

clmax = 1.8;
clmin = 0.2; 
cd0 = 0.03;    % with payload drag coefficient
cd0c = 0.025; % no external payload drag coefficient

% Function to calculate drag coefficient
% cd0 = dragBuildUp();

%Engine
TSFC = 13*1e-6;   %kg/(s.N)
prop_n = 0.8;


%mission parameters
reach_toc       = 2;  %hrs
cruise_alt      = 40000/3.281; %m
loiter_point    = 4000e3; %4000 km
v_cruise        = convvel(250,'kts','m/s');
v_loiter        = convvel(200,'kts','m/s');

target_roc      = cruise_alt/reach_toc/3600; %m/s
ld_climb        = 1/(2*sqrt(k*cd0)); %min l_d
cl_climb        = sqrt(pi*AR*cd0*e);

