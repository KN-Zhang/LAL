function [FP,FN] = plot_both_row(I1, I2, X, Y, VFCIndex, CorrectIndex, precise, recall,Fscore)
%   PLOT_MATCHES(I1, I2, X, Y, VFCINDEX, CORRECTINDEX)
%   considers correct indexes and indexes reserved by VFC, and then 
%   only plots the ture positive with blue lines, false positive with red
%   lines, false negative with green lines. For visibility, it plots at
%   most NUMPLOT (Default value is 50) matches proportionately.
%   
% Input:
%   I1, I2: Tow input images.
%
%   X, Y: Coordinates of intrest points of I1, I2 respectively.
%
%   VFCIndex: Indexes reserved by VFC.
%
%   CorrectIndex: Correct indexes.
%
%   See also:: VFC().

% Define the most matches to plot
NumPlot = 50;

n = size(X,1);
tmp=zeros(1, n);
tmp(VFCIndex) = 1;
tmp(CorrectIndex) = tmp(CorrectIndex)+1;
VFCCorrect = find(tmp == 2);
TruePos = VFCCorrect;   %Ture positive
tmp=zeros(1, n);
tmp(VFCIndex) = 1;
tmp(CorrectIndex) = tmp(CorrectIndex)-1;
FalsePos = find(tmp == 1); %False positive
tmp=zeros(1, n);
tmp(CorrectIndex) = 1;
tmp(VFCIndex) = tmp(VFCIndex)-1;
FalseNeg = find(tmp == 1); %False negative
tmp=zeros(1, n);
tmp(VFCIndex) = 1;
tmp(CorrectIndex) = tmp(CorrectIndex)-1;
TrueNeg=find(tmp == -1);

FP = FalsePos;
FN = FalseNeg;

all = 1:size(X,1);
TrueNeg = setdiff(all,FalsePos);
TrueNeg = setdiff(TrueNeg,TruePos);
TrueNeg = setdiff(TrueNeg,FalseNeg);

NumPos = length(TruePos)+length(FalsePos)+length(FalseNeg);
if NumPos > NumPlot
    t_p = length(TruePos)/NumPos;
    n1 = round(t_p*NumPlot);
    f_p = length(FalsePos)/NumPos;
    n2 = round(f_p*NumPlot);
    f_n = length(FalseNeg)/NumPos;
    n3 = round(f_n*NumPlot);
    t_n = length(TrueNeg)/NumPos;
    n4 = round(t_n*NumPlot);
else
    n1 = length(TruePos);
    n2 = length(FalsePos);
    n3 = length(FalseNeg);
    n4 = length(TrueNeg);
end

per = randperm(length(TruePos));
TruePos = TruePos(per(1:n1));
per = randperm(length(FalsePos));
FalsePos = FalsePos(per(1:n2));
per = randperm(length(FalseNeg));
FalseNeg = FalseNeg(per(1:n3));
per = randperm(length(TrueNeg));
TrueNeg = TrueNeg(per(1:n4));

% FalsePos = [FalsePos,24];

interval = 20;
WhiteInterval = 255*ones(size(I1,1), interval, size(I1,3));
figure
A=cat(2, I1, WhiteInterval, I2);
imshow(A) ;
hold on ;

% c=rand(size(X,1),3);
% for i=1:size(X,1)
% plot(X(i,1),X(i,2),'o','Color',c(i,:),'MarkerFaceColor',c(i,:),'MarkerSize',4);
% plot(Y(i,1)+size(I1,2)+interval,Y(i,2),'o','Color',c(i,:),'MarkerFaceColor',c(i,:),'MarkerSize',4);
% end

line([X(TruePos,1)'; Y(TruePos,1)'+size(I1,2)+interval], [X(TruePos,2)' ;  Y(TruePos,2)'],'linewidth', 1.2, 'color', 'b') ;%b
line([X(FalsePos,1)'; Y(FalsePos,1)'+size(I1,2)+interval], [X(FalsePos,2)' ;  Y(FalsePos,2)'],'linewidth', 1.2, 'color', 'r') ;%r
line([X(FalseNeg,1)'; Y(FalseNeg,1)'+size(I1,2)+interval], [X(FalseNeg,2)' ;  Y(FalseNeg,2)'],'linewidth', 1.2, 'color', 'g') ;%g
text(0,0,[sprintf('P: %.2f',precise) '  ' sprintf('R: %.2f',recall) '  ' sprintf('F: %.2f',Fscore)],'FontName','Times'...
     ,'Fontsize',35,'Color','w');

axis equal ;axis off  ; 
hold off
drawnow;



k = 0;
siz = size(I1);

figure


quiver(X(TrueNeg, 1), siz(1)-X(TrueNeg, 2), (Y(TrueNeg, 1)-X(TrueNeg, 1)), (-Y(TrueNeg, 2)+X(TrueNeg, 2)), k, 'k'), hold on%k
quiver(X(TruePos, 1), siz(1)-X(TruePos, 2), (Y(TruePos, 1)-X(TruePos, 1)), (-Y(TruePos, 2)+X(TruePos, 2)), k, 'b'), hold on%b
quiver(X(FalsePos, 1), siz(1)-X(FalsePos, 2), (Y(FalsePos, 1)-X(FalsePos, 1)), (-Y(FalsePos, 2)+X(FalsePos, 2)), k, 'r'), hold on%r
quiver(X(FalseNeg, 1), siz(1)-X(FalseNeg, 2), (Y(FalseNeg, 1)-X(FalseNeg, 1)), (-Y(FalseNeg, 2)+X(FalseNeg, 2)), k, 'g'), hold on%g



axis equal
axis([0 siz(2) 0 siz(1)]);
set(gca,'XTick',-2:1:-1)
set(gca,'YTick',-2:1:-1)






