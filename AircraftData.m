% Design 3 - Existing Data
% Author 460368355

clear
clc

[~,~,data] = xlsread('SimilarAircraftResearch');

columns = data(1,1:25);

uavOrNo = find(strcmp(columns,'UAV?'));
% Plot payload vs MTOW
iPayload = find(strcmp(columns,'Payload (kg)'));
iMTOW = find(strcmp(columns,'MTOW (kg)'));

for i = 2:33
    mtow = data{i,iMTOW};
    payload = data{i,iPayload};
    if ~(isnan(mtow) || isnan(payload))
        plot(mtow, payload, 'kx');
        hold on;
    end
end