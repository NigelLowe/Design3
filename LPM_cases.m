function [] = LPM_cases(X_Lon,X_Lat,T,xl,V)
% Aerodynamic Angles
alpha = rad2deg(X_Lon(2,:));
beta = rad2deg(X_Lat(1,:));
% Rotation rates
q = rad2deg(X_Lon(3,:));
p = rad2deg(X_Lat(2,:));
r = rad2deg(X_Lat(3,:));
% Attitude
theta = rad2deg(X_Lon(4,:));
phipsi = rad2deg(X_Lat(4:5,:));
% Velocities
u = X_Lon(1,:);
v = X_Lat(1,:).*V;
w = X_Lon(2,:);
% Plot response
set(groot,'defaultLineLineWidth',1.5) %Set the plotting line width to 2 thickness for all graph
%% Active states

subplot(1,2,1)
plot(T,theta);
xlabel('Time (s)');
ylabel('Pitch Angle (degrees)');
xlim(xl)
legend('cg1, 100knts', 'cg1, 180knts', 'cg2, 100knts', 'cg2, 180knts');
set(gca,'FontSize',15)
legend boxoff
grid on
grid minor
hold on

% Velocities
subplot(1,2,2)
plot(T,u);
xlabel('Time (s)');
ylabel('Velocity in x-body axis (m/s)');
xlim(xl)
set(gca,'FontSize',15)
legend('cg1, 100knts', 'cg1, 180knts', 'cg2, 100knts', 'cg2, 180knts');
legend boxoff
grid on
grid minor
hold on
end