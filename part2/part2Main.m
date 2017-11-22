%% Estimate fundamental matrix given two images and no matching points
%% Uses RANSAC
% Set flag for normalized and unset for unnormalized

clc;
close all;
clear;

flag=1; %% set 0 for unnormalized estimation

I1 = imread('library1.jpg');
I2 = imread('library2.jpg');

%% Estimation of Fundamental Matrix
[F,Inlier, matches] = fundamentalMatrixRANSAC(I1,I2,flag);
N=size(matches,1);
% Epipolar lines on image 2
L = (F * [matches(:,1:2) ones(N,1)]')';

% Finding epipolar lines closest to matches
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

% Finding endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

%% Visualization of epipolar lines

clf;
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

% Inliers
fprintf('Number of Inlier:'); disp(Inlier);

% Residuals
%A= dist2(matches(:,1:2),Th1(:,1:2));
%for i=1:1:size(matches,1)
 %   temp(i)=A(i,i);
%end
temp=sum(pt_line_dist)/length(pt_line_dist);
fprintf('Average Residual');disp(abs(temp));

