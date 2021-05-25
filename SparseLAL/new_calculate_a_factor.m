function neighbor_support_value=new_calculate_a_factor(X,Y,K_nn)%%V是场拟合出来的运动向量
%%
A=[X,X+Y,Y];gamma=10;
neighbor_support_index=FastCluster_J(A,K_nn,gamma);%%是个Nx3的矩阵
neighbor_support_index(find(neighbor_support_index(:,3)==0),3)=1;
N=size(X,1);
l_Y=sqrt(sum(Y.^2,2));%%是个列向量，记录了每一个运动向量的长度
cos_V=Y*Y'./(l_Y*l_Y');%%N*N的矩阵，角度的差别
K_nn=K_nn+1;
neighborX=knnsearch(X,X,'K',K_nn);%%idx的格式是N*K,且第一个是这个点本身
neighborX=neighborX(:,2:K_nn)';%%第j列记录了与第j个点最邻近的K个点的索引

neighborY=knnsearch(X+Y,X+Y,'K',K_nn);%%idx的格式是N*K
neighborY=neighborY(:,2:K_nn)';%%第j列记录了与第j个点最邻近的K个点的索引

neighborIndex = [neighborX; neighborY];
sorted_index= sort(neighborIndex);%按列升序
temp1 = diff(sorted_index);%%后一个减前一个
temp2 = (temp1 == zeros(size(temp1, 1), size(temp1, 2)));
temp2=[temp2;zeros(1,N)];
[rows_ini,cols]=ind2sub(size(temp2),find(temp2==1));
[sorted_temp_index,temp_index]=sort(cols);%%sorted_temp_index代表每个点的序号
rows_ini=rows_ini(temp_index);
rows=sorted_index(sub2ind(size(sorted_index),rows_ini,sorted_temp_index));
support_index=[sorted_temp_index rows];%%第一列是每个点的索引，第二列是这个点的支持点的索引

diff_length1=min([l_Y(support_index(:,1)),l_Y(support_index(:,2))],[],2)./max([l_Y(support_index(:,1)),l_Y(support_index(:,2))],[],2);
diff_cos1=cos_V(sub2ind(size(cos_V),support_index(:,1),support_index(:,2)));
diff_val1=diff_length1.*diff_cos1;

diff_length2=min([l_Y,l_Y(neighbor_support_index(:,3))],[],2)./max([l_Y,l_Y(neighbor_support_index(:,3))],[],2);
diff_cos2=cos_V(sub2ind(size(cos_V),(1:N)',neighbor_support_index(:,3)));
diff_val2=diff_length2.*diff_cos2;

neighbor_support_value1=accumarray(support_index(:,1),diff_val1,[N,1]);%%是个列向量
neighbor_support_value2=diff_val2;

neighbor_support_value1=neighbor_support_value1/(K_nn-1);
neighbor_support_value=0.7*neighbor_support_value1+0.3*neighbor_support_value2;

