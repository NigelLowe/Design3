clear all
clc
%% define parameters
% Aspect Ratio
AR = 10;
% gravity
g = 9.81;
% MTOW
MTOW = 22500;
% Air density
rho = 1.225;
% Cruise speed
V_cruise = 128.6;
% fuselage length
L_f = 16;
% Wing area
S_W = 48.4;
% Wing span
b = 22;
% Main wing
% Mean ac
wMac = 2.302;
% trailing edge angle = 13 (NASA 1015) thickness 0.152 - get location of ac
% in terms of mean aerodynamic chord
w_x_ac = 0.256*wMac;
% V-tail
% Mean ac
vMac = 1.89;
% trailing edge angle = 15 (NACA 2412) thickness 0.12 
t_x_ac = 0.248*vMac;
% Define horizontal lever arm (approximated)(Propeller in front of the
% fuselage 60%)
l_H = 6.72;
% Define vetical lever arm (approximated)(Propeller in front of the
% fuselage 60%)
l_V =  6.72;

% Deflection limits
% Elevator 20 deg down deflection, 25 deg up deflection
% Rudder 30 degrees left to right 


%% Sizing estimation
% Statistical horizontal tail volume coefficient 
C_H = 0.7;
% First estimate of horizontal tail area
% S_H = (C_H*S_W*wMac)/l_H;
% Final tail area
S_H = 0.2*S_W;


% Statistical vertical tail volume coefficient 
C_V = 0.065;
% First estimate of vertical tail area
S_V = (C_V*S_W*b)/l_V;

%% mass of horizontal tail plane
k_H = 1.1;
V_D = 1.4*V_cruise;
m_H = k_H*S_H*(62*((S_H^(0.2))*V_D)/1000*sqrt(cosd(4)) - 2.5);



%% Set up X_plot

%% Limit for Contorl Requirement 
% Coefficient of lift for main wing at low-speed control requirements
% Vstall at airfield, MTOW 22500 kg clmax 2(take-off rotation)
Vstall = 61;
V_R = 1.1*Vstall;
C_L = (MTOW*g)/((V_R^2)*(0.5*rho*S_W));
% x_cg percentage of fuselage
x_cg_bar = linspace(0.3,0.6,200);
% x_cg location in m
x_cg = L_f.*x_cg_bar;
% x_ac location in m
x_ac = 6.5  + w_x_ac;
% x_cg - x_ac
x_cgac = (x_cg - x_ac)/wMac;
% Assume for the moment that the z displacement of engine thrust is 0
C_ME = 0;
% Tail efficiency (will update to more accurate value)(literature value)
eff_H = 0.85;
% Coefficient of lift for horizontal tail
% Zero lift alpha for tail
% A_0H = -0.039;
% Downwash angle
% delta_1 = 2.28;
% delta_2 = 1.63;
% eps = (C_L/(pi*AR))*((1-(C_L/sqrt(C_L^2 + 1)))*delta_1 + (C_L/sqrt(C_L^2 + 1))*delta_2);
% C_L alpha (tail aerofoil naca 2412) (problem with xfoil
% prediction is the high reynolds number 1 e7)
C_Lah = 6.17;
% Angle of attack of 2 degrees
% alpha = deg2rad(15);
% i_H = deg2rad(-2.3);
% C_LH = C_Lah*(alpha + i_H - eps - A_0H);
% Tail (NACA0012) (need to be updated with more aerodynamic
% details)(literature value)
C_LH = -0.5;

% Moment about wing AC (wing aerofoil nasa 1015) (problem with xfoil
% prediction is the high reynolds number 1 e7)
C_MW = -0.107;
% Area ratio for take off rotation
Ratio_c = C_L/(C_LH*eff_H*(l_H/wMac)) * x_cgac + (C_MW + C_ME)/(C_LH*eff_H*(l_H/wMac));
% Area ratio for maximum CL
C_Lmax = 2.1;
% C_MWmax = (need an update once we have slats and flaps)
Ratio_f = C_Lmax/(C_LH*eff_H*(l_H/wMac)) * x_cgac + (C_MW + C_ME)/(C_LH*eff_H*(l_H/wMac));



%% Limit for Stability requirement
% C_L alpha (wing aerofoil nasa 1015) (problem with xfoil
% prediction is the high reynolds number 1 e7)
C_Law = 6.51;

% de_da (historical estimate 0.4 ~ 0.6)
de_da = 0.5;
% Area ratio for stability requirement (S_h/S_W)
Ratio_s = (C_Law*x_cgac)./((C_Lah*eff_H*(1-de_da)*((l_H/wMac)-x_cgac)));

% Neutral point limit
% V_bar = Ratio_s.*(l_H/wMac);
% Np = ((((x_ac)+eff_H*V_bar*(C_Law/C_Lah)*(1-de_da)))-x_ac)/wMac;
% C_Ma = 0.55;
Static_margin = 10/100; % Roskam II
a = C_Lah*eff_H*(1-de_da)*((l_H/wMac)-x_cgac);
% Area ratio for neutral point
Ratio_np = (Static_margin*C_Law-C_Law*x_cgac)./-a;

%% Tipback limit
% Choose lenght of landing gear from back of fuselage (Need to calculate
% how much weight is taken by the fron landing gear, ideal (20%)
% Movign rear landing gear foward means a longer landing gear which means a smaller
% V-tail height. This means a smaller v angle, which means a bigger
% horizontal tail area which means a very long tailplane span
B_L = 6.5;
 
% Tip back angle is calculated 10 deg, therefore, use 15 deg (literature suggested limiti) as limit
% from cg to strut 
% Tipback angle
tb_a = 15;
% Angle of attack requried for landing (15 deg) determines landing gear length (including
% wheel diameter)
LG_H = B_L*tand(tb_a); 


fprintf('Landing gear height: %.2f m\n',LG_H)

% Length between cg and rear landing gear
x = LG_H*tand(tb_a);
% fprintf('Length between cg and rear landing gear: %.2f m\n',x)
% tipback limit (where cg can be without the plane tipping back)
limit = (L_f - (B_L + x)- x_ac)/wMac;
% fprintf('cg location: %.2f m\n',L_f - (B_L + x))
% fprintf('cg location MAC: %.2f m\n',limit)
%%
% cg travel (0.4381 - 0.4692 at the moment) STANAG AMC.21 +/- 7% tolerance for cg travel 
cg_r1 = 0.4381;
cg_r2 = 0.4692;
cg1 = ((cg_r1.*L_f)-x_ac)/wMac;
cg2 = ((cg_r2.*L_f)-x_ac)/wMac;

%% limit ratio based on v tail angle
% Hangar height
H = 6;
% Fuselage height
F_H = 1.8;
% Allowable height in hangar considering landing gear height (not sure not
% high the tail is installed on the fuselage at the moment)
h_allow = 6 - LG_H;
fprintf('Allowable height for wing: %.2f m\n',h_allow)
% New V_tail angle
v = atand(S_V/S_H)    ; %(v angle of the tail)(dihedral angle)
fprintf('New V-tail angle: %.2f deg\n',v)
% Want angle between 20-30 degrees (too low would mean long span, too high
% would mean vertical tail height is too high)
v_angle = 40;
% S_H is now constraint by hangar height not x_cg
% S_H = S_V./tand(v_angle);
% New total tail area
New_total = S_H + S_V;
% New ratio
New_SH_SW = S_H./S_W;
fprintf('New S_h/S_w: %.2f \n',New_SH_SW)
fprintf('New horizontal wing area: %.2f m^2\n',S_H)
fprintf('New vertical wing area: %.2f m^2\n',S_V)
fprintf('New V_tail area: %.2f m^2\n',New_total)
% New v_tail span
New_b_H = (New_total*2)./(2.5+1.2);
fprintf('New V_tail span: %.2f m\n',New_b_H)
% Vertical tail height
VT_H = (New_b_H/2).*tand(v_angle);
fprintf('New V_tail height: %.2f m\n',VT_H)

%% Dimensioning for AVL
% root chord
cr = 3.03;
% tip chord
ct = 1.37;
% sweep angle
gamma = atand((cr - ct)/(New_b_H/2));
% Section 1
c1 = cr;
rv1 = 0.6;
lrv = 2.5-c1*rv1;
% Section 2
y2 = (New_b_H/2)*0.8;
x2 = y2*tand(gamma);
c2 = c1 - x2;
z2 = y2*tand(v_angle);
rv2 = lrv/c2;
% Section 3
y3 = (New_b_H)/2;
c3 = cr - ct;
z3 = y3*tand(v_angle);



%% Plot limits
% figure (1)
% plot(x_cgac,Ratio_c,'k')
% xlabel('C.G. location (% mac)');
% ylabel('S_h/S_w');
% ylim([0 0.8])
% xlim([-0.3 1.2])
% grid on
% grid minor
% hold on
% 
% plot(x_cgac,Ratio_f,'b')
% plot(x_cgac,Ratio_s,'r')
% plot(x_cgac,Ratio_np,'m')
% c = linspace(cg1,cg2,30);
% y = [];
% for i = 1:length(c)
%     y(i) = 0.17;
% end
% plot(c,y,'--','LineWidth',2)
% 
% 
% 
% 
% tb = xline(limit,'-',{'Tip back limit'});
% tb.LabelVerticalAlignment = 'bottom';
% tb.LabelHorizontalAlignment = 'left';
% tb.LineWidth = 2;
% tb.FontSize = 20;
% % Plot on X-plot to find solution
% hangar_limit = yline(0.467,'--',{'Hangar limit'});
% hangar_limit.LabelHorizontalAlignment = 'center';
% hangar_limit.LineWidth = 2;
% hangar_limit.FontSize = 20;
% hangar_limit.LabelVerticalAlignment = 'top';
% legend('Takeoff Rotation','Controllability at Cl_{max}','Stability Limit','Neutral Point','CG Travel');
% legend boxoff
% % 


