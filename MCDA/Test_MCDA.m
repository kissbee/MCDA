function  Test_MCDA(X, Y)
% Input:      X: data matrix
%             Y: true label
%             k: neighbor number
%             alpha,beta,gamma,theta1,theta2:parameters

% Output:     result: clustering result

% MSRC:
%    params:100000 5 0.0005 1 0.01
%    k:12

alpha=100000;
beta=5;
gamma=0.0005;
theta1=1;
theta2=0.01;

k = 12;

% V,n,c
V = size(X,2);      % view number
c = max(Y);         % class number
n = length(Y);      % object number

tic;

% Algorithem start.
F = MCDA(X,V,n,c,k,alpha,beta,gamma,theta1,theta2);

pre_label = kmeans(F,c,'maxiter',100,'replicates',10,'emptyaction','singleton');
result=Clustering8Measure(Y, pre_label);% result = [ACC nmi Purity Fscore Precision Recall AR Entropy];

% Display results
disp(result);

toc;
time=toc;
fprintf("run time：%f seconds",time);
