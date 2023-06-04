function  Test_MCDA(X, Y)
% Input:      X: data matrix
%             Y: true label
%             k: neighbor number
%             alpha,beta,gamma,theta1,theta2:parameters

% Output:     result: clustering result

%% k  * paras
% MSRC_v1	5	*50  0.5  5  1  0.5				
% MSRC	12	*100000  5  0.0005  1  0.01				
% HW	20	*100000  5000  1  0.00005  5  				
% NGs	7	*50 10000 0.1 0.05 0.0000005				
% WikipediaArticles	75	*100  100000  1  0.001  0.001 				
% Prokaryotic	33	*1000  10000  10  0.1  0.0001 				
% Uci_digit	14	*100000  1000  0.1  0.00005  0.001				
% Umist_mtv	5	*100  100  500000  1000  1000	
% Webkb	38	*1000  100  0.00001  1  0.1		

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
