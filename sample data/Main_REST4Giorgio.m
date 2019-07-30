% (r) REST scripts for UK_TMS_EEG_Exp
% Shiang Hu, 7/30/2019

clear; 
rRESTpath = 'D:\Andy\Desktop\scicol\Sample_of_EEG_file_with_xyz_channel_location\rREST_core';
eeglabpath = 'E:\OneDrive - Neuroinformatics Collaboratory\Scripting\Toolbox\eeglab';
addpath(eeglabpath); addpath(rRESTpath);
eeglab;  clc; clear; close all;

% prepare
eegname = 'Sample_of_EEG_file_with_xyz_channel_location.set';
eegpath = 'D:\Andy\Desktop\scicol\Sample_of_EEG_file_with_xyz_channel_location\';
lfname = 'Lead_Field_Subject_01_KB_23_04_2019_REAL_normal_-yxz.dat';

EEG = pop_loadset('filename',eegname,'filepath',eegpath);
data = double(EEG.data); % ref: FCz
data([29 30 64 65],:) = [];  % remove A1/2, V/HEOG
K = importdata(lfname)';
[Nc, Nt] = size(K);
H = rREST_Hsc(Nc,20); 
data = H*data;
EEG.nbchan = size(data,1);
EEG.chanlocs([29 30 64 65]) = [];

% REST
REST_ref = K*pinv(H*K, 0.05)*H;
EEG.data = REST_ref * data;
EEG.ref = 'REST';
EEG = pop_saveset(EEG,'filename',insertBefore(eegname, '.', '_REST_ref'),'filepath',eegpath);

% rREST
[data1, H1, L, s, lmd] = rREST_core(data,K); 
rREST_ref = pinv(L)*diag(1./(s+lmd))*H1';
EEG.data = rREST_ref * data1;
EEG.ref = 'rREST';
EEG = pop_saveset(EEG,'filename',insertBefore(eegname, '.', '_rREST_ref'),'filepath',eegpath);