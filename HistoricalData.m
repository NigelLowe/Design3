% Design 3 - Existing Data
% Author 460368355

clear
clc
close all

set(groot,'defaultLineLineWidth',2.0) %Set the plotting line width to 2 thickness for all graph
set(0,'defaultLineMarkerSize',30); % set default line marker size to 40
dataMax = 37;
[~,~,data] = xlsread('SimilarAircraftResearch');

columns = data(1,1:25);

uavOrNo = find(strcmp(columns,'UAV?'));


% Our aircraft
designPayload = 3500; % kg
designMtow = 37000; % kg
designEndurance = 46; % hr
designSpan = 73.5; % m
designRange = 0.9*(250*0.51444) * (designEndurance*3600) /1000; % km
designPower = 8202; % kW

% Payload vs MTOW - show higher mtow for higher payload
iPayload = find(strcmp(columns,'Payload (kg)'));
iMTOW = find(strcmp(columns,'MTOW (kg)'));
allMtow = [];
allPayload = [];
figure(1)
for i = 2:dataMax
    mtow = data{i,iMTOW};
    payload = data{i,iPayload};
    if ~(isnan(mtow) || isnan(payload))
        if strcmp(data{i,uavOrNo},'Yes')
            h1 = plot(mtow, payload, 'k.');
        else
            h2 = plot(mtow, payload, 'r.');
        end
        allMtow = [allMtow mtow];
        allPayload = [allPayload payload];
        hold on;
    end
end
p = polyfit(allMtow,allPayload,1);
x = linspace(0,max(allMtow),100);
y = p(1).*x + p(2);
plot(x,y,'k:');
h3 = plot(designMtow, designPayload, 'bx');
xlim([0 50000])
xlabel('MTOW (kg)');
ylabel('Payload (kg)');
set(gca,'FontSize',20)
set(gcf,'Position',[00 400 550 520])
legend([h1,h2,h3],'UAV','Manned','Our Solution', 'Location','northwest');
grid on


% Endurance vs MTOW - show lower mtow for good endurance. Opposite figure 1
iEndurance = find(strcmp(columns,'Flight Time (hr)'));
allMtow = [];
allEndurance = [];
figure(2)
for i = 2:dataMax
    mtow = data{i,iMTOW};
    endurance = data{i,iEndurance};
    if ~(isnan(mtow) || isnan(endurance))
        if strcmp(data{i,uavOrNo},'Yes')
            h1 = plot(mtow, endurance, 'k.');
        else
            h2 = plot(mtow, endurance, 'r.');
        end
        allMtow = [allMtow mtow];
        allEndurance = [allEndurance endurance];
        hold on;
    end
end
p = polyfit(allMtow,allEndurance,1);
x = linspace(0,max(allMtow),100);
y = p(1).*x + p(2);
plot(x,y,'k:');
h3 = plot(designMtow, designEndurance, 'bx');
xlim([0 50000])
xlabel('MTOW (kg)');
ylabel('Flight Time (hr)');
set(gca,'FontSize',20)
set(gcf,'Position',[00 400 550 520])
legend([h1,h2,h3],'UAV','Manned','Our Solution');
grid on


% Wing Span vs Range - show we want big span for good range
iSpan = find(strcmp(columns,'Wingspan (m)'));
iRange = find(strcmp(columns,'Max Range (km)'));
allRange = [];
allSpan = [];
figure(3)
for i = 2:dataMax
    span = data{i,iSpan};
    range = data{i,iRange};
    if ~(isnan(span) || isnan(range))
        if strcmp(data{i,uavOrNo},'Yes')
            h1 = plot(span, range, 'k.');
        else
            h2 = plot(span, range, 'r.');
        end
        allRange = [allRange range];
        allSpan = [allSpan span];
        hold on;
    end
end
p = polyfit(allSpan,allRange,1);
x = linspace(0,max(allSpan),100);
y = p(1).*x + p(2);
plot(x,y,'k:');
h3 = plot(designSpan, designRange, 'bx');
% ylim([0 15000]);
xlabel('Wingspan (m)');
ylabel('Range (km)');
set(gca,'FontSize',20)
set(gcf,'Position',[00 400 670 520])
legend([h1,h2,h3],'UAV','Manned','Our Solution', 'Location','southeast');
grid on

% uncomment this for comparing for cruise engine and propeller
% % Endurance vs Power - to show we want a small engine but have to
% % comprimise to takeoff (highlight C-2 Greyhound)
% iPower = find(strcmp(columns,'Power (kW)'));
% iModel = find(strcmp(columns,'Model'));
% LPPower = [];
% LPEndurance = [];
% figure(4)
% for i = 2:dataMax
%     endurance = data{i,iEndurance};
%     power = data{i,iPower}*1.34102;
%     if ~(isnan(endurance) || isnan(power))
%         if strcmp(data{i,uavOrNo},'Yes')
%             h1 = plot(endurance, power, 'k.');
%             if power < 1000
%                 LPPower = [LPPower power];
%                 LPEndurance = [LPEndurance endurance];
%             end
%         else
%             h2 = plot(endurance, power, 'r.');
%         end
%         hold on;
%     end
% end
% h3 = plot(designEndurance, designPower, 'bx');
% xlabel('Flight Time (hr)');
% ylabel('Horse Power (hp)');
% set(gca,'FontSize',20)
% set(gcf,'Position',[00 400 670 520])
% legend([h1,h2,h3],'UAV','Manned','Our Solution');
% grid on
% grid minor
% 
% 
% figure(5) % Zoom in of figure 4
% plot(LPEndurance, LPPower, 'k.')
% xlabel('Flight Time (hr)');
% ylabel('Horse Power (hp)');
% set(gca,'FontSize',20)
% set(gcf,'Position',[00 400 670 520])
% grid on
% grid minor


% Endurance vs Power - to show we want a small engine but have to
% comprimise to takeoff (highlight C-2 Greyhound)
iPower = find(strcmp(columns,'Power (kW)'));
iModel = find(strcmp(columns,'Model'));
allMtow = [];
allPower = [];
figure(6)
for i = 38:47
    mtow = data{i,iMTOW};
    power = data{i,iPower};
    if ~(isnan(mtow) || isnan(power))
        if strcmp(data{i,uavOrNo},'Helicopter')
            h1 = plot(mtow, power, 'k.');
        end
        allMtow = [allMtow mtow];
        allPower = [allPower power];
        hold on;
    end
end
p = polyfit(allMtow,allPower,1);
x = linspace(0,max(allMtow),100);
y = p(1).*x + p(2);
plot(x,y,'k:');
h2 = plot(designMtow, designPower, 'bx');
xlabel('MTOW (kg)');
ylabel('Power (kW)');
set(gca,'FontSize',20)
set(gcf,'Position',[00 400 670 520])
legend([h1,h2],'Helicopter/VTOL','Our Solution', 'Location','northwest');
grid on

