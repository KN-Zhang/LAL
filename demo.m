clear;close all;
startup; %% setup vlfeat
%% loading data
addpath('./data/');
I1=imread('nc_0837.jpg');  %% Candidate frame
I2=imread('nc_1961.jpg');  %% Query image
load nc_0837vs1961.mat;   %% The putative set constructed in advance. If not exist, users could create it based on the following codes.

%%%% no initial correspondences
% SiftThreshold = 1.5; % no smaller than 1
% [X, Y] = sift_match(I1, I2, SiftThreshold);

%% LAL & LAL*
MismatchRemoval_method='LAL';
addpath([MismatchRemoval_method '/']);
[nX, nY, normal]=norm2(X,Y);
if ~exist('conf', 'var'), conf = []; end
conf = VFC_init(conf);
K_nn=6;

tic;           
switch MismatchRemoval_method
    case 'LAL'
        VecFld=LAL(nX, nY-nX, K_nn, conf);
    case 'SparseLAL'
        VecFld=SparseLAL(nX, nY-nX, K_nn, conf);
end
index=VecFld.VFCIndex;
toc
rmpath([MismatchRemoval_method '/']);
%% Evaluation
if ~exist('CorrectIndex', 'var'), CorrectIndex = index; end
[precise, recall, corrRate] = evaluate(CorrectIndex, index, size(X,1));
Fscore=2*precise*recall/(precise+recall)

%% Plot results 
plot_both_row(I1, I2, X, Y, index, CorrectIndex, precise, recall, Fscore);
rmpath('./data/');