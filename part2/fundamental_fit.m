%% Function to estimate fundamental matrix given the matching points
function [F] = fundamental_fit(matches,flag)
%% Input
%  F-> Fundamental Matrix
%  matches-> matched inliers/points->[4*4]
%  flag-> Set flag for normalization, unset for estimation without normalization

%% Output
%  F-> Output Fundamental Matrix 

a = matches(:,1:2); %% Points from Image 1
b = matches(:,3:4); %% Points from Image 2

%% Normalization
if (flag == 1)
    % for first image
    a(:,3)=1;
    mat=a;
    MEAN = mean(mat');
    MAXIMUM = max(mat');
    T1 = [inv(MAXIMUM(1)), 0, -MEAN(1)*inv(MAXIMUM(1));0, inv(MAXIMUM(2)), -MEAN(2)*inv(MAXIMUM(2)); 0, 0, 1];
    match = T1*mat';
    a=match';
    % for second image
    b(:,3)=1;
    mat=b;
    MEAN = mean(mat');
    MAXIMUM = max(mat');
    T2 = [inv(MAXIMUM(1)), 0, -MEAN(1)*inv(MAXIMUM(1));0, inv(MAXIMUM(2)), -MEAN(2)*inv(MAXIMUM(2)); 0, 0, 1];
    match = T2*mat';
    b=match';
end

%% A Matrix Calculation and Fundamental matrix generation

x=a(:,1);
y=a(:,2);
X=b(:,1);
Y=b(:,2);
A=[];

for i=1:1:length(x)
    A=[A;X(i)*x(i),X(i)*y(i),X(i),Y(i)*x(i),Y(i)*y(i),Y(i),x(i),y(i), 1];
    %A=[A;X(i)*x(i),x(i)*Y(i),x(i),y(i)*X(i),Y(i)*y(i),y(i),X(i),Y(i), 1];
    %%Another variant of A formula
end

[U,S,V]=svd(A);
V=V(:,end);
V=reshape(V,[3,3]);
[A,B,C]=svd(V,0);
B(end,end)=0; %% Setting Rank=2 by setting the minimum element as zero
F=A*B*C';

if (flag==1)
    F=T2'*F*T1;
end


end
