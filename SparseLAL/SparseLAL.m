function VecFld = SparseLAL(X, Y, K_nn, conf)
% VFC  Vector Field Consensus
%   VECFLD = SPARSEVFC(X, Y, CONF)
%   learns vector field from random samples with outliers. It is a fast
%   version of VFC. Even fast then FastVFC(), and has better performance
%   than FastVFC.
%   
% Input:
%   X, Y: Original data.
%
%   conf: configuration for VFC. See VFC_init().
%
% Output:
%   VecFld: A structure type value which contains X, Y, beta, V, C, P,
%       VFCIndex. Where V = f(X), P is the posterior probability and 
%       VFCIndex is the indexes of inliers which found by VFC.
%
%   See also:: VFC_init(), VFC(), FastVFC().

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

fprintf('Start mismatch removal...\n');
[N, D]=size(Y);

% Calculate Support Value
neighbor_support_value=new_calculate_a_factor(X,Y,K_nn);
an=zeros(size(X,1),1);
sn=neighbor_support_value;
an(sn>0)=-sn(sn>0)+1;%%max(sn)
an(sn<=0)=1;
bn=(1-neighbor_support_value)./(conf.a*neighbor_support_value);
bn(bn<=0)=1;

% Construct kernel matrix K
M=10;
tmp_X = unique(X, 'rows'); idx = randperm(size(tmp_X,1)); 

idx = idx(1:min(M,size(tmp_X,1)));ctrl_pts=tmp_X(idx,:);%%ctrl_pts相当于新选出来的基底
K=con_K(ctrl_pts,ctrl_pts,conf.beta);
U = con_K(X, ctrl_pts, conf.beta);
M = size(ctrl_pts,1);%%选出来的基底个数

% Initialization
V=zeros(N,D); iter=1;  tecr=1; C=zeros(M,D); E=1; 
sigma2 = sum(sum((Y-V).^2))/(N*D); gamma = conf.gamma;
%%
while (iter < conf.MaxIter) && (tecr > conf.ecr) && (sigma2 > 1e-8) 
    % E-step. 
    E_old=E;
    [P, E]=get_P( Y, V, an, bn, sigma2, gamma);
    E=E+conf.lambda/2*trace(C'*K*C);
    tecr=abs((E-E_old)/E);
%     fprintf('iterate: %dth, gamma: %f, the energy change rate: %f, sigma2=%f\n', iter, gamma, tecr, sigma2);

    % M-step. Solve linear system for C.
    P = max(P, conf.minP);
    C=(U'.*repmat(P', [M, 1])*U+conf.lambda*sigma2*K)\(U'.*repmat(P', [M, 1])*Y);

    % Update V and sigma^2
    V=U*C;
    Sp=sum(P);
    sigma2=sum(P'*sum((Y-V).^2, 2))/(Sp*D);

    % Update gamma
    numcorr = length(find(P > conf.theta));
    gamma=numcorr/size(X, 1);
    if gamma > 0.95, gamma = 0.95; end
    if gamma < 0.05, gamma = 0.05; end
    
    iter=iter+1;
end

%%
VecFld.X = ctrl_pts;
VecFld.Y = Y;
VecFld.beta = conf.beta;
VecFld.V=V;
VecFld.C=C;
VecFld.P = P;
VecFld.VFCIndex = find(P > conf.theta);
VecFld.ninliers=numel(VecFld.VFCIndex);
VecFld.sigma2 = sigma2;

disp('Removing outliers succesfully completed.');


%%%%%%%%%%%%%%%%%%%%%%%%
function K=con_K(x,y,beta)
% CON_K constructs the kernel K, 
%   where K(i, j) = k(x, y) = exp(-beta*||x-y||^2).

[n, d]=size(x); [m, d]=size(y);

K=repmat(x,[1 1 m])-permute(repmat(y,[1 1 n]),[3 2 1]);
K=squeeze(sum(K.^2,2));
K=-beta * K;
K=exp(K);
