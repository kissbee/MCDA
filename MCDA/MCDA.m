function F  = MCDA(X,V,n,c,k,alpha,beta,gamma,theta1,theta2)
%% =========================Initialization=========================

% Lv1 Lv2
for iv = 1:V
    S = constructS_PNG(X{iv}, k, 1); %%  input: dv * n
    S1= (S+S')/2;
    D1 = diag(sum(S1));
    Lv1{iv} = D1-S1; % first-order Laplacian matrix
    
    S2=zeros(n,n);
    for i=1:n
        for j=1:n
            S2(i,j) = S1(:,i)'*S1(:,j);
        end
    end
    D2 = diag(sum(S2));
    Lv2{iv} = D2-S2; % secend-order Laplacian matrix
end

% Lv Fv Rv
for iv = 1:V
    Lv{iv} = theta1*Lv1{iv}+theta2*Lv2{iv};
    
    [Fviv, ~, ~]=eig1(Lv{iv}, c, 0);
    Fv{iv} = Fviv;
    
    Rv{iv} = ones(c, c);
end

% Wv dv pv qv mv
for iv=1:V
    dv(iv) = size(X{iv}, 1); % dimension of each view
    Wv{iv} = ones(dv(iv), c);
    pv(iv) = 1.0/V;
    qv(iv) = 1.0/V;
    mv(iv) = 1.0/V;
end

% L  F
L = ones(n, n);
F = ones(n, c);
F = orth(F);%%%

maxIter = 100;

%% =========================Normalization=========================
% Normalization
for i = 1:V
    fea = X{i}';
    fea = mapminmax(fea,0,1);  % n * dv
    X{i}= fea'; % dv * n
end

%% =========================Optimization=========================
iter = 1;
while iter<=maxIter
    
    % Updating pv, qv, mv
    for iv = 1:V
        zv(iv) = norm(X{iv}'* Wv{iv}-F,'fro');
    end
    for iv = 1:V
        qv(iv) = 0.5/norm(Lv{iv}-L,'fro');
        mv(iv) = 0.5/norm(Fv{iv}*Rv{iv}-F,'fro');
        pv(iv) = zv(iv)/sum(zv);
    end
    
    % Updating Wv
    for iv = 1:V
        I = eye(dv(iv));
        Wv{iv} = 1.0/pv(iv) * inv(1.0/pv(iv)*X{iv}*X{iv}'+alpha*I)*X{iv}*F;
    end
    
    % Updating L
    up = zeros(n,n);
    down = 0;
    for iv = 1:V
        up = up + qv(iv)*Lv{iv};
        down = down + qv(iv);
    end
    L = (2.0*beta*up-F*F')/(2.0*beta*down);
    
    % Updating F
    [~,y] = eig(L);
    m = diag(y);
    lamdaMax = max(m); % Find the maximum characteristic of L
    
    C = zeros(n, c);
    C_ = zeros(n, c);
    for iv = 1:V
        C = C + 1.0/pv(iv)*X{iv}'*Wv{iv};
        C_ = C_ + mv(iv)*Fv{iv}*Rv{iv};
    end
    C = C+gamma*C_;
    
    F1 = ones(n,c);
    I=eye(n);
    while 1
        M  = 2*(lamdaMax*I-L)*F+2*C;
        [U_,~,V_] = svds(M);
        F = U_ * V_';
        if(norm(F-F1))<1e-3 % convergence
            break;
        end
        F1 = F;
    end
    
    % Updating  Rv
    for iv=1:V
        Rv{iv} = inv(Fv{iv}'*Fv{iv})*Fv{iv}'*F;
    end
    
    % Calculating obj value
    tempObj = 0;  % Item 1
    for iv = 1:V
        tempObj = tempObj + 1.0 / pv(iv) * norm(X{iv}' * Wv{iv} - F, 'fro')^2;
    end
    
    sumWv = 0;  % Item 2
    for iv = 1:V
        sumWv = sumWv + norm(Wv{iv},'fro')^2;
    end
    tempObj = tempObj + alpha * sumWv;
    
    sumLL = 0;  % Item 3
    for iv = 1:V
        sumLL = sumLL + qv(iv) * norm(Lv{iv} - L, 'fro')^2;
    end
    tempObj = tempObj + beta * sumLL;
    
    sumFF = 0;  % Item 4
    for iv = 1:V
        sumFF = sumFF + mv(iv) * norm(Fv{iv}*Rv{iv} - F, 'fro')^2;
    end
    tempObj = tempObj + gamma * sumFF;
    
    tempObj = tempObj + trace(F' * L * F);  % Item 5
    
    obj(iter) = tempObj;  % record 
    
    
    % Convergence checking
%     if iter > 1  && abs(obj(iter)-obj(iter-1))/obj(iter-1) < 1e-3
%         break;
%     end
    
    iter = iter+1;
end

%% =============================Plot=============================

% Plot convergence curve
plot(obj);
xlabel('The number of iterations');
ylabel('Object function value');

end