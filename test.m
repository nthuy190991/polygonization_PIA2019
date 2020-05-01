clear
clc
close all

load('example/P_U_shape')
% load('example/P_L_shape')
% load('example/P_rect')

GP=[]; % GP empty for the original algorithm
Vspec=1;
Mspec=1;
resampling=0; % without resampling
verbose=0; 

[vertices,~,S,Psub]=polygonization_dutter(P,GP,Vspec,Mspec,resampling,verbose);

figure
plot(P(:,1),P(:,2),'b.')
hold on, grid on, axis equal
plot(vertices(1,:),vertices(2,:),'r')
