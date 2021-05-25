function [precision, recall, corrRate] = evaluate(CorrectIndex, VFCIndex, siz)
%   [PRECISION, RECALL, CORRRATE] = EVALUATE(CORRECTINDEX, VFCINDEX, SIZ)
%   evaluates the performence of VFC with precision and recall.
%
% Input:
%   CorrectIndex, VFCIndex: Ground truth indexes and indexes preserved by VFC.
%
%   siz: Number of initial matches.
%
% Output:
%   precision, recall, corrRate: Precision and recall of VFC, percentage of
%       initial correct matches.
%
%   See also:: VFC().

% Authors: Jiayi Ma (jyma2010@gmail.com)
% Date:    04/17/2012

if length(VFCIndex)==0
    VFCIndex = 1:siz;
end

VFCCorrect = intersect(VFCIndex, CorrectIndex);
NumCorrectIndex = length(CorrectIndex);
NumVFCIndex = length(VFCIndex);
NumVFCCorrect = length(VFCCorrect);

corrRate = NumCorrectIndex/siz;
precision = NumVFCCorrect/NumVFCIndex;
recall = NumVFCCorrect/NumCorrectIndex;

fprintf('\ncorrect correspondence rate in the original data: %d/%d = %f\n', NumCorrectIndex, siz, corrRate);
fprintf('precision rate: %d/%d = %f\n', NumVFCCorrect, NumVFCIndex, precision);
fprintf('recall rate: %d/%d = %f\n', NumVFCCorrect, NumCorrectIndex, recall);