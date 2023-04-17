clc;
clear memory;
warning('off')
addpath('data')
addpath('func')

dataname = 'MSRC_v1';         % set dataset name here. 
disp(['--test data:'  dataname ])

load(dataname);            % load data

Test_MCDA(X,Y);            % test dataset dimension: dv * n.
