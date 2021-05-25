%function [P,P_diag,sigma2_new, E]=get_P(X,K,Y,V, sigma2 ,gamma, a)
function [P, E]=get_P(Y,V,an,bn, sigma2 ,gamma)
% GET_P estimates the posterior probability and part of the energy.
%%
% neighbor_support_value=new_calculate_a_factor(X,Y,K);%%������
% op_neighbor_support_value=max(neighbor_support_value)-neighbor_support_value;
% an=zeros(size(X,1),1);
% sn=neighbor_support_value;
% an(sn>0)=-sn(sn>0)+max(sn);
% an(sn<=0)=max(sn);
% an=mapminmax(an',0,1);
% an=an';

D = size(Y, 2);
temp1 = exp(-an.*sum((Y-V).^2,2)/(2*sigma2));%%������
temp2 = (2*pi*sigma2)^(D/2)*(1-gamma)./(gamma);%%һ����


% bn=(max(neighbor_support_value)-neighbor_support_value)./(a*neighbor_support_value);

% P=(neighbor_support_value.*temp1)./(neighbor_support_value.*temp1+op_neighbor_support_value.*temp2);%%�˴���P�Ǹ�����
P=(temp1)./(temp1+bn.*temp2);%%�˴���P�Ǹ�����

% E=P'*sum((Y-V).^2,2)/(2*sigma2)+sum(P)*log(sigma2)*D/2;%%(3.16)ǰ���е��෴��
E=sum(P.*an.*sum((Y-V).^2,2))/(2*sigma2)+sum(P)*log(sigma2)*D/2;
%%
% mixture_pai=calculate_a_factor(X,Y,K);%%NxN����
% N=size(Y,1);
% D = size(Y, 2);
% YY=permute(repmat(Y,[1,1,N]),[3 2 1]);
% VV=repmat(V,[1,1,N]);
% tmp=squeeze(sum((YY-VV).^2,2));%%��һ���ǵ�һ����������N����֮��ľ����ƽ��
% tmp=tmp';%%��һ���ǵ�һ����������N����֮��ľ����ƽ��
% tmp1=mixture_pai.*exp(-tmp/(2*sigma2));
% tmp2 = (2*pi*sigma2)^(D/2)*(1-gamma)./(gamma*a);
% P=tmp1./repmat(sum(tmp1,2)+tmp2,[1 N 1]);%%�˴���P�Ǹ�����
% P_diag=sum(P,2);
% MP=sum(P(:));
% sigma2_new=trace(P*tmp')/(MP*D);
% E=trace(P*tmp')/(2*sigma2)+ MP*log(sigma2)*D/2-MP*log(1-gamma)-(N-MP)*log(gamma);%%����Ҫ��С������������
