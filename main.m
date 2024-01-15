clc; clear; close all;

%% load data
load('T.mat');

%% 7SD
[Ps_7SD,Pd_7SD,Pv_7SD,Ph_7SD,Pod_7SD,Pcd_7SD,Pmd_7SD] = sevenSD(T);

% % 保存数据，用ENVI画图
% Save_data_to_ENVI(Ps_7SD,1,'float','./output/Ps_7SD');
% Save_data_to_ENVI(Pd_7SD,1,'float','./output/Pd_7SD');
% Save_data_to_ENVI(Pv_7SD,1,'float','./output/Pv_7SD');
% Save_data_to_ENVI(Ph_7SD,1,'float','./output/Ph_7SD');
% Save_data_to_ENVI(Pod_7SD,1,'float','./output/Pod_7SD');
% Save_data_to_ENVI(Pcd_7SD,1,'float','./output/Pcd_7SD');
% Save_data_to_ENVI(Pmd_7SD,1,'float','./output/Pmd_7SD');



