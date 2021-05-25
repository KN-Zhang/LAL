%function [P,P_diag,sigma2_new, E]=get_P(X,K,Y,V, sigma2 ,gamma, a)
function [P, E]=get_P(Y,V,an,bn, sigma2 ,gamma)
% GET_P estimates the posterior probability and part of the energy.
%%
% neighbor_support_value=new_calculate_a_factor(X,Y,K);%%列向量
% op_neighbor_support_value=max(neighbor_support_value)-neighbor_support_value;
% an=zeros(size(X,1),1);
% sn=neighbor_support_value;
% an(sn>0)=-sn(sn>0)+max(sn);
% an(sn<=0)=max(sn);
% an=mapminmax(an',0,1);
% an=an';

D = size(Y, 2);
temp1 = exp(-an.*sum((Y-V).^2,2)/(2*sigma2));%%列向量
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)./(gamma);%%一个数


% bn=(max(neighbor_support_value)-neighbor_support_value)./(a*neighbor_support_value);

% P=(neighbor_support_value.*temp1)./(neighbor_support_value.*temp1+op_neighbor_support_value.*temp2);%%此处的P是个向量
P=(temp1)./(temp1+bn.*temp2);%%此处的P是个向量

% E=P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;%%(3.16)前两行的相反数
E=sum(P.*an.*sum((Y-V).^2,2))/(2*sigma2)+sum(P)*log(sigma2)*D/2;
%%
% mixture_pai=calculate_a_factor(X,Y,K);%%NxN矩阵
% N=size(Y,1);
% D = size(Y, 2);
% YY=permute(repmat(Y,[1,1,N]),[3 2 1]);
% VV=repmat(V,[1,1,N]);
% tmp=squeeze(sum((YY-VV).^2,2));%%第一列是第一个点与其他N个点之间的距离的平方
% tmp=tmp';%%第一行是第一个点与其他N个点之间的距离的平方
% tmp1=mixture_pai.*exp(-tmp/(2*sigma2));
% tmp2 = (2*pi*sigma2)^(D/2)*(1-gamma)./(gamma*a);
% P=tmp1./repmat(sum(tmp1,2)+tmp2,[1 N 1]);%%此处的P是个矩阵
% P_diag=sum(P,2);
% MP=sum(P(:));
% sigma2_new=trace(P*tmp')/(MP*D);
% E=trace(P*tmp')/(2*sigma2)+ MP*log(sigma2)*D/2-MP*log(1-gamma)-(N-MP)*log(gamma);%%这是要最小化的能量函数
