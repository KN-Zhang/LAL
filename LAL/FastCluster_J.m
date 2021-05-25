function neighbor_support_index=FastCluster_J(A,K_nn,gamma)%%���ص���������
% clear all  
% close all  
% addpath('./data/');
% % load('motion_fields_77.mat');
% % N_points=size(motion_fields_77.X,1);
% load('malaga_0001vs0523.mat');
N_points=size(A,1);
% disp('The only input needed is a distance matrix file')  
% disp('The format of this file should be: ')  
% disp('Column 1: id of element i')  
% disp('Column 2: id of element j')  
% disp('Column 3: dist(i,j)')  
  
%% ���ļ��ж�ȡ����  
% mdist=input('name of the distance matrix file (with single quotes)?\n');  
xx=zeros(N_points*N_points,3);
tmp=[1:N_points]';
xx(:,1)=reshape((repmat(tmp,1,N_points,1))',N_points*N_points,1);%111111222222333333....
xx(:,2)=repmat(tmp,N_points,1,1);%1234123412341234....
D=myDist(A,gamma);%%��ʱ��D�Ǹ�NxN�ľ���,A��Nx6�ľ���
% D=pdist2(A,A);
% D=D';
xx(:,3)=D(:);

% D=pdist(motion_fields_77.V);
% D=squareform(D);
% xx(:,3)=D(:);
% disp('Reading input distance matrix')  
% xx=load(mdist);  


% ND=max(xx(:,2));  
% NL=max(xx(:,1));  
% if (NL>ND)  
%   ND=NL;  %% ȷ�� DN ȡΪ��һ�������ֵ�еĽϴ��ߣ���������Ϊ���ݵ�����  
% end  
ND=N_points;
NL=N_points;
N=size(xx,1); %% xx ��һ��ά�ȵĳ��ȣ��൱���ļ�����������������ܸ�����  
  
%% ��ʼ��Ϊ��  
% for i=1:ND  
%   for j=1:ND  
%     dist(i,j)=0;  
%   end  
% end  
 
%% ���� xx Ϊ dist ���鸳ֵ��ע������ֻ���� 0.5*DN(DN-1) ��ֵ�����ｫ�䲹����������  
%% ���ﲻ���ǶԽ���Ԫ��  
% for i=1:N  
%   ii=xx(i,1);  
%   jj=xx(i,2);  
%   dist(ii,jj)=xx(i,3);  
%   dist(jj,ii)=xx(i,3);  
% end  
  dist=D;
%% ȷ�� dc  
  
percent=5.0;  
% fprintf('average percentage of neighbours (hard coded): %5.6f\n', percent);  
  
position=round(N*percent/100); %% round ��һ���������뺯��  
sda=sort(xx(:,3)); %% �����о���ֵ����������  
dc=sda(position);  
  
%% ����ֲ��ܶ� rho (���� Gaussian ��)  
  
% fprintf('Computing Rho with gaussian kernel of radius: %12.6f\n', dc);  
  
%% ��ÿ�����ݵ�� rho ֵ��ʼ��Ϊ��  
% for i=1:ND  
%   rho(i)=0.;  
% end  
  rho=zeros(ND,1);
  dist00 =exp(-(dist./dc).^2);
% Gaussian kernel  
for i=1:ND-1  
  for j=i+1:ND  
     rho(i)=rho(i)+dist00(i,j);  
     rho(j)=rho(j)+dist00(i,j);  
  end  
end  
  
% "Cut off" kernel  
%for i=1:ND-1  
%  for j=i+1:ND  
%    if (dist(i,j)<dc)  
%       rho(i)=rho(i)+1.;  
%       rho(j)=rho(j)+1.;  
%    end  
%  end  
%end  
  
%% ������������ֵ���������ֵ�����õ����о���ֵ�е����ֵ  
maxd=max(max(dist));   
  
%% �� rho ���������У�ordrho ������  rhoΪ��ľֲ��ܶ�
[rho_sorted,ordrho]=sort(rho,'descend');  
   
%% ���� rho ֵ�������ݵ�  
delta(ordrho(1))=-1.;  
nneigh(ordrho(1))=0;  
  
%% ���� delta �� nneigh ����  
for ii=2:ND  
   delta(ordrho(ii))=maxd;  
   for jj=1:ii-1  
     if(dist(ordrho(ii),ordrho(jj))<delta(ordrho(ii)))  
        delta(ordrho(ii))=dist(ordrho(ii),ordrho(jj));  
        nneigh(ordrho(ii))=ordrho(jj);   
        %% ��¼ rho ֵ��������ݵ����� ordrho(ii) ��������ĵ�ı�� ordrho(jj)  
     end  
   end  
end  
  
%% ���� rho ֵ������ݵ�� delta ֵ  
delta(ordrho(1))=max(delta(:));  
rho_delta=rho.*delta;  
% ����ͼ  
  
% disp('Generated file:DECISION GRAPH')   
% disp('column 1:Density')  
% disp('column 2:Delta')  
  
% fid = fopen('DECISION_GRAPH', 'w');  
% for i=1:ND  
%    fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));  
% end  
  
%% ѡ��һ��Χס�����ĵľ���  
% disp('Select a rectangle enclosing cluster centers')  
  
%% ÿ̨�����������ĸ�����ֻ��һ����������Ļ�����ľ������ 0  
%% >> scrsz = get(0,'ScreenSize')  
%% scrsz =  
%%            1           1        1280         800  
%% 1280 �� 800 ���������õļ�����ķֱ��ʣ�scrsz(4) ���� 800��scrsz(3) ���� 1280  
% scrsz = get(0,'ScreenSize');  
  
%% ��Ϊָ��һ��λ�ã��о���û����ô auto �� :-)  
% figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);  
  
%% ind �� gamma �ں��沢û���õ�  
% for i=1:ND  
% %   ind(i)=i;   
%   rho_delta(i)=rho(i)*delta(i);  
% end  

[~,sorted_rho_delta_index]=sort(rho_delta,'descend');



%% ���� rho �� delta ����һ����ν�ġ�����ͼ��  
  
% subplot(2,1,1)  
% tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');  
% title ('Decision Graph','FontSize',15.0)  
% xlabel ('\rho')  
% ylabel ('\delta')  
  
% subplot(2,1,1)  
% rect = getrect(figure(1));   
%% getrect ��ͼ��������ȡһ���������� rect �д�ŵ���  
%% �������½ǵ����� (x,y) �Լ����ؾ��εĿ�Ⱥ͸߶�  
% rhomin=rect(1);  
% deltamin=rect(2); %% ���߳������Ǹ� error������ 4 ��Ϊ 2 ��!  

  
%% ��ʼ�� cluster ����  
NCLUST=0;  
  
%% cl Ϊ������־���飬cl(i)=j ��ʾ�� i �����ݵ�����ڵ� j �� cluster  
%% ��ͳһ�� cl ��ʼ��Ϊ -1  
% for i=1:ND  
%   cl(i)=-1;  
% end  
  cl=zeros(ND,1)-1;
%% �ھ���������ͳ�����ݵ㣨���������ģ��ĸ���  icl��¼�˾������ĺ�
% for i=1:ND  
%   if ( (rho(i)>rhomin) && (delta(i)>deltamin))  
%      NCLUST=NCLUST+1;  
%      cl(i)=NCLUST; %% �� i �����ݵ����ڵ� NCLUST �� cluster  
%      icl(NCLUST)=i;%% ��ӳ��,�� NCLUST �� cluster ������Ϊ�� i �����ݵ�  
%   end  
% end  
NCLUST=round(ND/(2*K_nn));
% NCLUST=5;
cl(sorted_rho_delta_index(1:NCLUST))=1:NCLUST; 
for i=1:NCLUST
   icl(i)= sorted_rho_delta_index(i);
end

for i=1:NCLUST
   [~,sorted_D_index]=sort(D(icl(i),icl));
   ncl(i)=icl(sorted_D_index(2));   
end
% icl_ncl=[icl;ncl];
% fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);  
  
% disp('Performing assignation')  
  
%% ���������ݵ���� (assignation)  
for i=1:ND  
  if (cl(ordrho(i))==-1)  
    cl(ordrho(i))=cl(nneigh(ordrho(i)));  
  end  
end  
%% �����ǰ��� rho ֵ�Ӵ�С��˳�����,ѭ��������, cl Ӧ�ö��������ֵ��.   
 
%% ������ε㣬halo��δ���Ӧ���Ƶ� if (NCLUST>1) ��ȥ�ȽϺð�  
% for i=1:ND  
%   halo(i)=cl(i);  
% end  
  
% if (NCLUST>1)  
  
  % ��ʼ������ bord_rho Ϊ 0,ÿ�� cluster ����һ�� bord_rho ֵ  
%   for i=1:NCLUST  
%     bord_rho(i)=0.;  
%   end  
  
  % ��ȡÿһ�� cluster ��ƽ���ܶȵ�һ���� bord_rho  
%   for i=1:ND-1  
%     for j=i+1:ND  
%       %% �����㹻С��������ͬһ�� cluster �� i �� j  
%       if ((cl(i)~=cl(j))&& (dist(i,j)<=dc))  
%         rho_aver=(rho(i)+rho(j))/2.; %% ȡ i,j �����ƽ���ֲ��ܶ�  
%         if (rho_aver>bord_rho(cl(i)))   
%           bord_rho(cl(i))=rho_aver;  
%         end  
%         if (rho_aver>bord_rho(cl(j)))   
%           bord_rho(cl(j))=rho_aver;  
%         end  
%       end  
%     end  
%   end  
  
  %% halo ֵΪ 0 ��ʾΪ outlier  
%   for i=1:ND  
%     if (rho(i)<bord_rho(cl(i)))  
%       halo(i)=0;  
%     end  
%   end  
  
% end  
  
%% ��һ����ÿ�� cluster  
% for i=1:NCLUST  
%   nc=0; %% �����ۼƵ�ǰ cluster �����ݵ�ĸ���  
%   nh=0; %% �����ۼƵ�ǰ cluster �к������ݵ�ĸ���  
%   for j=1:ND  
%     if (cl(j)==i)   
%       nc=nc+1;  
%     end  
%     if (halo(j)==i)   
%       nh=nh+1;  
%     end  
%   end  
% %   fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);  
% end  
  
% hot=colormap;  
% for i=1:NCLUST  
%    ic=int8((i*64.)/(NCLUST*1.));  
%    subplot(2,1,1)  
%    hold on  
%    plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',hot(ic,:),'MarkerEdgeColor',hot(ic,:));  
% end  
% subplot(2,1,2)  
% disp('Performing 2D nonclassical multidimensional scaling')  
% Y1 = mdscale(dist, 2, 'criterion','metricstress');  
% plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');  
% title ('2D Nonclassical multidimensional scaling','FontSize',15.0)  
% xlabel ('X')  
% ylabel ('Y')  
% for i=1:ND  
%  A(i,1)=0.;  
%  A(i,2)=0.;  
% end  
% for i=1:NCLUST  
%   nn=0;  
%   ic=int8((i*64.)/(NCLUST*1.));  
%   for j=1:ND  
%     if (halo(j)==i)  
%       nn=nn+1;  
%       A(nn,1)=Y1(j,1);  
%       A(nn,2)=Y1(j,2);  
%     end  
%   end  
%   hold on  
%   plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',hot(ic,:),'MarkerEdgeColor',hot(ic,:));  
% end
% nclncl=zeros(1,ND);
% nclncl=ncl(cl);
% nclncl(find(halo==0))=0;
neighbor_support_index=[(1:ND)',cl,ncl(cl)'];
%for i=1:ND  
%   if (halo(i)>0)  
%      ic=int8((halo(i)*64.)/(NCLUST*1.));  
%      hold on  
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));  
%   end  
%end  
% faa = fopen('CLUSTER_ASSIGNATION', 'w');  
% disp('Generated file:CLUSTER_ASSIGNATION')  
% disp('column 1:element id')  
% disp('column 2:cluster assignation without halo control')  
% disp('column 3:cluster assignation with halo control')  
% for i=1:ND  
%    fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));  
% end