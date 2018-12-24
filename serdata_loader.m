% serdata_loader
% just load a bunch of serdata

disp('loading data...wait...')
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_drug0.mat')
serdata_drug0 = serdata;
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_ps0.mat')
serdata_ps0 = serdata;
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_sc1.mat')
serdata_sc1 = serdata;
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_sc0.mat')
serdata_sc0 = serdata;
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_ps1.mat')
serdata_ps1 = serdata;
load('Z:\Katsuhisa\serotonin_project\LFP_project\Data\serdata_drug1.mat')
serdata_drug1 = serdata;
disp('loaded!')
clearvars serdata