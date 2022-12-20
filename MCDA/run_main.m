clc;
clear memory;
warning('off')
addpath('data')
addpath('func')

dataname = 'MSRC';         % set dataset name here. 
disp(['--test data:'  dataname ])

load(dataname);            % load data

for i = 1:size(X,2)
    X{i} = X{i}';
end

Test_MCDA(X,Y);            % test dataset dimension: dv * n.
