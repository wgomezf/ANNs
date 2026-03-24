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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [X,FEATS] = SHAPE_FEATS(BW) computes 12 morphological features from the 
%   binary blob of a breast tumor BW. X is a numeric vector with the feature 
%   values and FEATS is a cell vector with the name of the features in the 
%   same order as in X.
function [x,feats] = shape_feats(BW)
    % Padding to ensure tumor pixels has background pixels
    BW2 = padarray(BW,[1 1],0,'both');
    %---------------------------------------------------------------------
    % Normalized Residual Value
    BW_props = regionprops(BW2,'ConvexHull','Perimeter','Centroid','Area');
    CH = roipoly(BW2,BW_props.ConvexHull(:,1),BW_props.ConvexHull(:,2));
    NRV = bwarea(xor(BW2,CH))/bwarea(CH); 
    %---------------------------------------------------------------------
    % Compute the tumor's equivalent ellipse
    Pbw = regionprops(BW,'Area','Centroid','Perimeter');
    A = Pbw.Area;
    xc = Pbw.Centroid(1); yc = Pbw.Centroid(2);
    [y,x] = find(BW);
    [xx,yy] = meshgrid(1:size(BW,2),1:size(BW,1));
    % Calculate the second-order moments of the original binary object
    Sxx = (1/A)*sum((x-xc).^2);
    Syy = (1/A)*sum((y-yc).^2);
    Sxy = (1/A)*sum((x-xc).*(y-yc));
    % Calculate the coefficients of the general equation of the ellipse
    coef = (1/(4*Sxx*Syy - Sxy^2))*[Syy -Sxy;-Sxy Sxx];
    a = coef(1,1); b = coef(1,2); c = coef(2,2);
    % Calculate the equivalent ellipse of the binary object
    E = (a*(xx-xc).^2 + 2*b*(xx-xc).*(yy-yc) + c*(yy-yc).^2) < 1;
    % Calculate the perimeter, and minor and major axes of the equivalent ellipse
    E_props = regionprops(E,'Perimeter','MinorAxisLength','MajorAxisLength');
    %---------------------------------------------------------------------
    % Overlap ratio
    OR  = bwarea(BW&E)/bwarea(BW|E);	
    %---------------------------------------------------------------------
    % Elliptic-normalized circumference
    ENC = E_props.Perimeter/Pbw.Perimeter; 
    S = bwmorph(BW,'skel','inf');
    S = bwareaopen(S,5);
    %---------------------------------------------------------------------
    % Elliptic-normalized skeleton
    ENS = sum(S(:))/E_props.Perimeter; 
    %---------------------------------------------------------------------
    %Long axis to short axis ratio
    LS = E_props.MajorAxisLength/E_props.MinorAxisLength; 
    Emax = E_props.MajorAxisLength; % Major axis
    Emin = E_props.MinorAxisLength; % Minor axis
    %---------------------------------------------------------------------
    % Compactness
    R = 1 - ((4*pi*A)/(Pbw.Perimeter^2)); 
    %---------------------------------------------------------------------
    % Shape class
    % Parameterization of tumor and ellipse contours
    junk = bwboundaries(BW);
    cBW  = junk{1};
    yBW  = cBW(:,1); xBW = cBW(:,2);
    junk = bwboundaries(E);
    cE  = junk{1};
    yE  = cE(:,1); xE = cE(:,2);
    % Unit vectors of the tumor contour
    rBW = [xBW-xc yBW-yc];
    nBW = sqrt(sum(rBW.^2,2));
    uBW = rBW./(repmat(nBW,1,2)+eps);
    % Unit vectors of the ellipse contour
    rE = [xE-xc yE-yc];
    nE = sqrt(sum(rE.^2,2));
    uE = rE./(repmat(nE,1,2)+eps);
    % Distance between unit vectors
    D1 = dist(uBW,uE');
    [~,ind] = min(D1,[],2);
    % Correspondence between tumor points and ellipse points 
    % with the closest orientation
    mdE = cE(ind,:);
    % Euclidean distance
    D2 = sqrt((cBW(:,1)-mdE(:,1)).^2+(cBW(:,2)-mdE(:,2)).^2);
    SC = sum(D2)/Pbw.Perimeter;
    %---------------------------------------------------------------------
    % Proportional distance
    S1 = bwperim(BW);
    S2 = bwperim(E);
    avdis1 = averagedist(S1,S2);
    avdis2 = averagedist(S2,S1);
    PD = 100*((avdis1 + avdis2)/(2*sqrt(bwarea(E)/pi)));
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Tumor's major axis orientation
    theta = 0.5*atan((2*Sxy)/(Sxx-Syy));
    if (Sxx > Syy) && (theta > 0)
        theta = 1*theta;
    elseif (Sxx > Syy) && (theta < 0)
        theta = -1*theta;
    elseif (Sxx < Syy) && (theta > 0)
        theta = pi/2 - theta;
    elseif (Sxx < Syy) && (theta < 0)
        theta = pi/2 - (-1*theta);
    end
    O = theta*180/pi;
    %++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    % Depth to width ratio
    [yBW,xBW] = find(BW);
    xBWmax = max(xBW); xBWmin = min(xBW);
    yBWmax = max(yBW); yBWmin = min(yBW);
    DWR = (yBWmax-yBWmin)/(xBWmax-xBWmin);
    %---------------------------------------------------------------------
    % Feature vector
    x = [NRV OR ENC ENS LS Emax Emin R SC PD O DWR];
    feats = {'NRV','OR','ENC','ENS','LS','AX_MX','AX_MN','ROUND','SC','PD','Angle','DWR'};
end
%---------------------------------------------------------------------
% Average distance between two contours
function avdis = averagedist(cs,cr)
    [lseg,cseg] = find(cs);
    [lreal,creal] = find(cr);
    [Lseg,Lreal] = meshgrid(lseg,lreal);
    [Cseg,Creal] = meshgrid(cseg,creal);
    dist = sqrt((Lseg-Lreal).^2+(Cseg-Creal).^2);
    clear Lseg Lreal Cseg Creal
    d = min(dist);
    avdis = mean(d);
end