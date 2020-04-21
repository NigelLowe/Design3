function [A_lat,B_lat] = LatMatrix_Trim_folda(q_bar,V,m)
%  Vortex Lattice Output -- Total Forces
% 
g = 9.81;
Ixx = 0.222e6;
Iyy = 0.8549e5;
Izz = 0.2933e6;
Ixz = 657.2;
%  Configuration: Albatross                                                   
%      # Surfaces =   4
%      # Strips   =  80
%      # Vortices = 800
% 
%   Sref =  48.400       Cref =  2.3100       Bref =  22.000    
%   Xref =  7.1983       Yref = 0.11397E-07   Zref =  3.6695    
  S =  48.400;       
  c_bar =  2.3055 ;      
  b =  22.000  ;  
%  Standard axis orientation,  X fwd, Z down         
% 
%  Run case: Trim                                    
% 
  Alpha =   deg2rad(-0.27766 )    ;
% pb/2V =  -0.00000     p'b/2V =  -0.00000
  Beta  =   deg2rad(0.00000)    ;
% qc/2V =   0.00000
%   Mach  =     0.000     rb/2V =  -0.00000     r'b/2V =   0.00000
% 
%   CXtot =  -0.03017     Cltot =  -0.00000     Cl'tot =  -0.00000
%   CYtot =  -0.00000     Cmtot =   0.00000
%   CZtot =  -1.00035     Cntot =   0.00000     Cn'tot =   0.00000
% 
%   CLtot =   1.00000
%   CDtot =   0.04006
%   CDvis =   0.00000     CDind = 0.0400646
%   CLff  =   0.99304     CDff  = 0.0391243    | Trefftz
%   CYff  =  -0.00000         e =    0.8023    | Plane  

   df            =   deg2rad(0.00000);
   da         =   deg2rad(0.00000);
   de        =  deg2rad(-7.98598);
   dr          =  deg2rad(-0.00046);

%  ---------------------------------------------------------------

%  Stability-axis derivatives...

%                              alpha                beta
%                   ----------------    ----------------
%  z' force CL |    
 CLa =   6.208350;    
 CL = CLa*Alpha;
 CLb =  -0.000000;
%  y  force CY |    
 CYa =   0.000000;    
 CYb =  -0.201900;
%  x' mom.  Cl'|    
 Cla =  -0.000000;    
 Clb =  -0.087071;
%  y  mom.  Cm |    
 Cma =  -1.539911;    
 Cmb =   0.000000;
%  z' mom.  Cn'|    
 Cna =  -0.000003;    
 Cnb =   0.056166;

%                      roll rate  p'      pitch rate  q'        yaw rate  r'
%                   ----------------    ----------------    ----------------
%  z' force CL |    
 CLp =   0.000000;    
 CLq =  20.390684;    
 CLr =  -0.000001;
%  y  force CY |   
 CYp =  -0.004384;   
 CYq = -0.013148;    
 CYr =   0.178713;
% x' mom.  Cl'|    
Clp =  -0.543304;   
Clq =   0.000000;    
Clr =   0.273330;
%  y  mom.  Cm |    
 Cmp =   0.000003;    
 Cmq = -23.480890;    
 Cmr =   0.000001;
%  z' mom.  Cn'|    
 Cnp =  -0.084253;    
 Cnq =  -0.000004;    
 Cnr =  -0.058550;

%                   Flap         df     Aileron      da     Elevator     de     Rudder       dr 
%                   ----------------    ----------------    ----------------    ----------------
% z' force CL |   
CLdf =   0.044698;   
CLda =  -0.000000;   
CLde =   0.011758;   
CLdr =   0.017637;
%  y  force CY |   
 CYdf =   0.000000;   
 CYda =  -0.000861;   
 CYde =  -0.000000;   
 CYdr =   0.001769;
%  x' mom.  Cl'|   
 Cldf =   0.000000;   
 Clda =  -0.003293;   
 Clde =   0.000000;   
 Cldr =   0.000011;
%  y  mom.  Cm |   
 Cmdf =   0.002714;
 Cmda =   0.000000;   
 Cmde =  -0.036463;   
 Cmdr =  -0.054695;
%  z' mom.  Cn'|   
 Cndf =  -0.000000;   
 Cnda =   0.000196 ;  
 Cnde =   0.000000 ; 
 Cndr =  -0.000635;
 %  Trefftz drag| 
 Cddf =  0.002986 ; 
 Cdda =   0.000000; 
 Cdde =  -0.000321; 
 Cddr =  -0.000481;
%  span eff.   |    edf =   0.010572    eda =  -0.000000    ede =   0.022171    edr =   0.033256



%  Neutral point  Xnp =   7.711332

%  Clb Cnr / Clr Cnb  =   0.338261    (  > 1 if spirally stable )
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
A_lat = [ A11 A12 A13 A14 0 ;
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
B_lat = [ B11 B12 ;
B21 B22 ;
B31 B32 ;
0 0 ;
0 0 ] ;
end
