function VecFld = LAL(X, Y, K_nn, conf)
% VFC  Vector Field Consensus
%   VECFLD = VFC(X, Y, CONF)
%   learns vector field from random samples with outliers.
%   
% Input:
%   X, Y: Original data.
%
%   conf: Configuration for VFC. See VFC_init().
%
% Output:
%   VecFld: A structure type value which contains X, Y, beta, V, C, P,
%       VFCIndex. Where V = f(X), P is the posterior probability and 
%       VFCIndex is the indexes of inliers which found by VFC.
%

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

fprintf('Start mismatch removal...\n');
[N, D]=size(Y); 

% Calculate Support Value
neighbor_support_value=new_calculate_a_factor(X,Y,K_nn);
an=zeros(size(X,1),1);
sn=neighbor_support_value;
an(sn>0)=-sn(sn>0)+1;
an(sn<=0)=1;
bn=(1-neighbor_support_value)./(conf.a*neighbor_support_value);
bn(bn<=0)=1;
  
% Construct kernel matrix K
K=con_K(X,X,conf.beta);

% Initialization
V = zeros(N,D); iter = 1;  tecr = 1; C = zeros(N,D); E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D); gamma = conf.gamma;
%%
while (iter < conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) 
    % E-step. 
      E_old=E;
      [P, E]=get_P( Y, V, an, bn, sigma2, gamma);
      E = E + conf.lambda/2*trace(C'*K*C);
      tecr = abs((E-E_old)/E);

    % M-step. Solve linear system for C.
      P = max(P, conf.minP);
      inv_P = diag(1./(P));
      C = (K+conf.lambda*sigma2*inv_P)\Y;

    % Update V and sigma^2
      V = K*C;
      Sp = sum(P);
      sigma2 = sum(an.*P.*sum((Y-V).^2, 2))/(Sp*D);  %%P is a column vector

    % Update gamma
      numcorr = length(find(P > conf.theta));
      gamma = numcorr/size(X, 1);
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter = iter+1;
end

%%
% Fix the value of gamma, redo the EM process. 
V = zeros(N,D); iter = 1;  tecr = 1; C = zeros(N,D); E = 1; 
sigma2 = sum(sum((Y-V).^2))/(N*D);
while (iter < conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) 
    % E-step.
    E_old = E;
    [P, E]=get_P( Y, V, an, bn, sigma2, gamma);  
    E = E+conf.lambda/2*trace(C'*K*C);
    tecr = abs((E-E_old)/E);

    % M-step. Solve linear system for C.
    P = max(P, conf.minP);
    inv_P = diag(1./(P));
    C = (K+conf.lambda*sigma2*inv_P)\Y;

    % Update V and sigma^2
    V = K*C;
    Sp = sum(P);
    sigma2 = sum(an.*P.*sum((Y-V).^2, 2))/(Sp*D);   %%P is a column vector

    iter = iter+1;
end
%%
VecFld.X = X;
VecFld.Y = Y;
VecFld.beta = conf.beta;
VecFld.V = V;
VecFld.C = C;
VecFld.P = P;
VecFld.VFCIndex = find(P > conf.theta);
VecFld.ninliers=numel(VecFld.VFCIndex);

disp('Outlier removal succesfully completes.');


%%%%%%%%%%%%%%%%%%%%%%%%
function K=con_K(x,y,beta)
% CON_K constructs the kernel K, 
%   where K(i, j) = k(x, y) = exp(-beta*||x-y||^2).
[n, d]=size(x); [m, d]=size(y);
K=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
K=squeeze(sum(K.^2,2));
K=-beta * K;
K=exp(K);

