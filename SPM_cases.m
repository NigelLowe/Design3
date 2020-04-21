function [] = SPM_cases(X_Lon,X_Lat,T,xl,V)
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
%% Active States
% figure(3)
subplot(3,1,1)
plot(T,alpha);
xlabel('Time (s)');
ylabel('\alpha (degrees)');
xlim(xl)
legend('cg1, 100knts', 'cg1, 180knts', 'cg2, 100knts', 'cg2, 180knts');
legend('Location', 'eastoutside')
set(gca,'FontSize',13)
legend boxoff
grid on
grid minor
hold on

subplot(3,1,2)
% figure(4)
plot(T,theta);
xlabel('Time (s)');
ylabel('\theta (degrees)');
xlim(xl)
set(gca,'FontSize',13)
legend('cg1, 100knts', 'cg1, 180knts', 'cg2, 100knts', 'cg2, 180knts');
legend('Location', 'eastoutside')
legend boxoff
grid on
grid minor
hold on

subplot(3,1,3)
% figure (5)
plot(T,q);
xlabel('Time (s)');
ylabel('q (deg/s)');
xlim(xl)
set(gca,'FontSize',13)
legend('cg1, 100knts', 'cg1, 180knts', 'cg2, 100knts', 'cg2, 180knts');
legend('Location', 'eastoutside')
legend boxoff
grid on
grid minor
hold on

end
