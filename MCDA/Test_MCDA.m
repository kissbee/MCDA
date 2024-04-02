function  Test_MCDA(X, Y)
% Input:      X: data matrix
%             Y: true label
%             k: neighbor number
%             alpha,beta,gamma,theta1,theta2:parameters

% Output:     result: clustering result

k = 5;

alpha=50;
beta=0.5;
gamma=5;
theta1=1;
theta2=0.5;

% V,n,c
V = size(X,2);      % view number
c = max(Y);         % class number
n = length(Y);      % object number

tic;

%% Algorithem start.
F = MCDA(X,V,n,c,k,alpha,beta,gamma,theta1,theta2);

for i=1:20
    pre_label = kmeans(F,c,'maxiter',100,'replicates',10,'emptyaction','singleton');
    result(i,:)=Clustering8Measure(Y, pre_label);% result = [ACC nmi Purity Fscore Precision Recall AR Entropy];
end
res = [mean(result);std(result)];

%% Display results
disp(res);

toc;
time=toc;
fprintf("run time：%f seconds",time);
