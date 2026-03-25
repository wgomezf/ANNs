%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: Wilfrido Gómez-Flores (Cinvestav)     %
% e-mail: wgomez@cinvestav.mx                   %
% Date:   february 2026                         %
% Subject: Compute morphological features       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;

% Opt: 1 - Malignant tumor and 0 - Benign tumor
opt = 1; 
if opt
    fileId = '119xzr8vpJP4JI2lwTPa-SQJhpIbiiER0';
else
    fileId = '1cAAMpyeRoRw11w_7ZZeAI4iEuBrLA7jf';
end
fileName = 'bus.mat';
fileUrl = sprintf('https://drive.google.com/uc?export=download&id=%s', fileId);
websave(fileName, fileUrl);
load(fileName)

% Display images
figure('color','w')
subplot 121;
imshow(I);
title(sprintf('BUS image with %s tumor',info.tumor.pathology));
subplot 122;
imshow(BW);
title('Tumor mask');

% Calculate shape features
[x,feats] = shape_feats(BW);
for i = 1:numel(x)
    fprintf('x[%d] - %s: %0.2f\n',i, feats{i},x(i));
end