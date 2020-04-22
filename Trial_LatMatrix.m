function [A_Lat,B_Lat] = Trial_LatMatrix(q_bar,V,Cd,m,AR)

%  ---------------------------------------------------------------
g = 9.81;
Ixx = 0.222e6;
Iyy = 0.8549e5;
Izz = 0.2933e6;
Ixz = 657.2;
% m = 22500;   
% AR = 10;
% V = 128;
% rho = 0.39;
% q_bar = 0.5*rho*V^2;
% Cd = 0.034;





%  Vortex Lattice Output -- Total Forces
%
%  Configuration: Albatross
%      # Surfaces =   4
%      # Strips   =  80
%      # Vortices = 800

file = 'Trim_fold2';
[Geo, Aero, Con] = readAVLForces(file);
        
e =    0.7675;     
k = 1/(pi*AR*e);


        S = Geo.S;   % Platform Area (m^2)
    
    
        c_bar = Geo.c ;   % Chord Length (m)
    
     
        b = Geo.b ;   % Wing Span (m)
    
    
    % Aerodynamic Data (Reference CG:  % mac)
    %find Alpha
     
        Alpha = Aero.Alpha;
    
    %find CL
      
        CL = Aero.CL;
    
    %find CD
      
        CD = Aero.CD ;
    
    %find Mach number
     
        MN = Aero.MN ;
    
    
    % Lt Coefficients
     
        CLa = Aero.CLa ;
    
   
        CLq = Aero.CLq ;
    
     
        CLde = Aero.CLde ;
    
      
        CLr = Aero.CLr ;
    
    
    % Side Force Coefficients

        CYb = Aero.CYb ;
    
 
        CYp = Aero.CYp ;
    
     
        CYr = Aero.CYr ;
    

        CYda = Aero.CYda ;
    

        CYdr = Aero.CYdr;
    
    
    % M Moment Coefficients
    %find CM

        CM = Aero.CM;
    

        Cma = Aero.Cma ;
    
  
        Cmq = Aero.Cmq ;
    

        Cmde = Aero.Cmde;
    
    
    % N Moment Coefficients
     
        Cnb = Aero.Cnb ;
    

        Cnp = Aero.Cnp ;
    

        Cnr = Aero.Cnr ;
    

        Cnda = Aero.Cnda;
    

        Cndr = Aero.Cndr;
    
    
    % L Moment Coefficients

        Clb = Aero.Clb;
    
 
        Clp = Aero.Clp ;
    

        Clr = Aero.Clr;
    

        Clda = Aero.Clda;
  

        Cldr = Aero.Cldr;
    
    
    % Control surfaces

        df = Con.df;
    
   
        da = Con.da;
    
 
        de = Con.de;
    

        dr = Con.dr;
        
        Cddf = Con.df;
        Cdda = Con.da;
        Cdde = Con.de;
        Cddr = Con.dr;
    

% Forces in X direction
Xu = -(2/(m*V))*q_bar*S*(Cd);
Cda = 0.007; % Wing airfoil property
Xa = (1/m)*q_bar*S*(CL + CLde*de - Cda);
% Xwdot = 0 ;
% Xq = 0;
Xde = -(1/m)*q_bar*S*Cdde  ;

% Forces in Z direction
Zu = -(2/(m*V))*q_bar*S*CL ;
Za = -(1/m)*q_bar*S*(CLa + Cd + Cdde*de);

de_da = 1.5*k*CLa;
l_t = 6.72;
CZadot = (Cmq*de_da)/(l_t/c_bar);
% Zwdot = q_bar*S*c_bar*CLadot /(2*m*V^2) ;
% Zq = q_bar*S*c_bar*CLq /(2*m*V) ;
Zde = -(1/m)*q_bar*S*CLde  ;

% M moments
CMu = 0;
Mu = (1/Iyy)*q_bar*S*c_bar*CMu;
Ma = (1/Iyy)*q_bar*S*c_bar*Cma ;
Cmadot = Cmq*de_da;
Madot = q_bar*S*(c_bar^2)*Cmadot / (2*Iyy*V);
Mq = q_bar*S*(c_bar^2)*Cmq /(2*Iyy*V) ;
Mde = (1/Ixx)*q_bar*S*c_bar*Cmde ;


% Inertia Ratios
% A1 = Ixz / Ixx ;
% B1 = Ixz / Izz ;

% Steady level Flight Condition
theta1 = 0 ;
%% Matricies
% Side forces 
Yb = (1/m)*q_bar*S*CYb ;
Yp = q_bar*S*b*CYp /(2*m*V) ;
Yr = q_bar*S*b*CYr /(2*m*V) ;
Yda = (1/m)*q_bar*S*CYda  ;
Ydr = (1/m)*q_bar*S*CYdr;

% L moments
Lb = (1/Ixx)*q_bar*S*b*Clb;
Lp = q_bar*S*(b^2)*Clp / (2*Ixx*V) ;
Lr = q_bar*S*(b^2)*Clr/ (2*Ixx*V) ;
Lda = (1/Ixx)*q_bar*S*b*Clda ;
Ldr = (1/Ixx)*q_bar*S*b*Cldr  ;

% N moments
Nb = (1/Izz)*q_bar*S*b*Cnb;
Np = q_bar*S*(b^2)*Cnp / (2*Izz*V) ;
Nr = q_bar*S*(b^2)*Cnr / (2*Izz*V) ;
Nda = (1/Izz)*q_bar*S*b*Cnda ;
Ndr = (1/Izz)*q_bar*S*b*Cndr;
% Inertia Ratios
A1 = Ixz / Ixx ;
B1 = Ixz / Izz ;


% Steady level Flight Condition
theta1 = 0 ;

%% Matricies
%  States Matrix (beta form) (Assume yaw moment derivative with respect to
%  sideslip due to thrust is zero)
% u = V cos(theta) where theta was assumed to be close to 0 at level flight
% which allowed small angle approximation 
A11 = Yb/(V) ;
A12 = Yp/(V) ;
A13 = (Yr/(V)) -1;
A14 = (g*cos(theta1)) / (V) ;
A21 = (Lb+(A1*Nb))/(1 - A1*B1 ) ;
A22 = (Lp+A1*Np)/(1 - A1*B1 ) ;
A23 = ( Lr+A1*Nr )/(1 - A1*B1 ) ;
A31 = (Nb+B1*Lb)/(1 - A1*B1 ) ;
A32 = (Np+B1*Lp)/(1 - A1*B1 ) ;
A33 = (Nr+B1*Lr )/(1 - A1*B1 ) ;
A_Lat = [ A11 A12 A13 A14 0 ;
          A21 A22 A23 0 0 ;
          A31 A32 A33 0 0 ;
          0 1 tan(theta1) 0 0 ;
          0 0 sec(theta1) 0 0 ] ;
      
% Control Matrix
B11 = Yda/V;
B12 = Ydr/V;
B21 = (Lda+A1*Nda)/(1 - A1*B1 ) ;
B22 = (Ldr+A1*Ndr )/(1 - A1*B1 ) ;
B31 = (Nda+B1*Lda )/(1 - A1*B1 ) ;
B32 = (Ndr+B1*Ndr )/(1 - A1*B1 ) ;
B_Lat = [ B11 B12 ;
B21 B22 ;
B31 B32 ;
0 0 ;
0 0 ] ;
end
% 






