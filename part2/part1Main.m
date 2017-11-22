%% Estimate fundamental matrix given matching points in 2 images - normalized and unnormalized
% Set flag for normalized and unset for unnormalized

clc;
close all;
clear;

flag=0; %% set 0 for unnormalized estimation

I1 = imread('house1.jpg');
I2 = imread('house2.jpg');
matches = load('house_matches.txt'); 
P1 = load('house1_camera.txt');
P2 = load('house2_camera.txt');
N = size(matches,1);

%imshow([I1 I2]); hold on;
%plot(matches(:,1), matches(:,2), '+r');
%plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
%line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');
%pause;

%% Estimation of Fundamental Matrix
F = fundamental_fit(matches,flag);

% Epipolar lines on image 2
L = (F * [matches(:,1:2) ones(N,1)]')';

% Finding epipolar lines closest to matches
L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); % rescale the line
pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);

if (flag==2)
    coor=find(abs(pt_line_dist)<15);
    matches=matches(coor,:);
    N = size(matches,1);
    L = (F * [matches(:,1:2) ones(N,1)]')';
    L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3);
    pt_line_dist = sum(L .* [matches(:,3:4) ones(N,1)],2);
    closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);
    
    
end
Average_Residue=(sum(pt_line_dist))/N;

% Finding endpoints of segment on epipolar line (for display purposes)
pt1 = closest_pt - [L(:,2) -L(:,1)] * 10; % offset from the closest point is 10 pixels
pt2 = closest_pt + [L(:,2) -L(:,1)] * 10;

%% Visualization of epipolar lines

clf;
imshow(I2); hold on;
plot(matches(:,3), matches(:,4), '+r');
line([matches(:,3) closest_pt(:,1)]', [matches(:,4) closest_pt(:,2)]', 'Color', 'r');
line([pt1(:,1) pt2(:,1)]', [pt1(:,2) pt2(:,2)]', 'Color', 'g');

fprintf('Average Residue:');disp(abs(Average_Residue));
