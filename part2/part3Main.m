%% This function estimates and visualizes the camera center and trianguation

clc;
close all;
clear;

% Read inputs
I1 = imread('house1.jpg');
I2 = imread('house2.jpg');
matches = load('house_matches.txt'); 
P1 = load('house1_camera.txt');
P2 = load('house2_camera.txt');
N = size(matches,1);

figure(1);
imshow([I1 I2]); hold on;
plot(matches(:,1), matches(:,2), '+r');
plot(matches(:,3)+size(I1,2), matches(:,4), '+r');
line([matches(:,1) matches(:,3) + size(I1,2)]', matches(:,[2 4])', 'Color', 'r');

%% Triangulation & Camera center

[Val,C1,C2] = triangulation1(matches,P1,P2); % this is a function that you should write
%L = (F * [matches(:,1:2) ones(N,1)]')'; % transform points from 
% the first image to get epipolar lines in the second image

Th1=[P1*Val']';
Th2=[P2*Val']';

%% Visualization of camera center and triangulation
figure(3)
axis equal;
scatter3(Val(:,1),Val(:,2),Val(:,3),'red');
hold on;
scatter3(C1(1),C1(2),C1(3),'*blue');
scatter3(C2(1),C2(2),C2(3),'ogreen');
title('Triangulation and Camera Center');

%% Residual calculation



%% 1st Image
A= dist2(matches(:,1:2),Th1(:,1:2));
for i=1:1:size(matches,1)
    temp=A(i,i);
end
temp=mean(temp);
fprintf('Residual of 1st image');disp(temp);
%% 2nd Image
A= dist2(matches(:,1:2),Th1(:,1:2));
for i=1:1:size(matches,1)
    temp(i)=A(i,i);
end
temp=sum(temp)/length(temp);
fprintf('Average Residual of 1st image');disp(temp);

B= dist2(matches(:,1:2),Th2(:,1:2));
for i=1:1:size(matches,1)
    temp(i)=B(i,i);
end
temp=sum(temp)/length(temp);
fprintf('Average Residual of 2nd image');disp(temp);


