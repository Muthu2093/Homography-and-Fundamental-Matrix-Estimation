%% Triangulation and Camera Center Estimation
function [Val,C1,C2]=triangulation1(matches,Point1,Point2)
%% Inputs
% matches-> Matching points
% Point1-> Camera co-ordinates of image 1
% Point2-> Camera co-ordinates of image 2

%% Outputs
% Val -> Triangulation
% C1 -> Camera center of image 1
% C2 -> Camera center of image 2

%I1=imread('/Users/muthuvel/Documents/MATLAB/hw3/data/part2/house1.jpg');
%I2=imread('/Users/muthuvel/Documents/MATLAB/hw3/data/part2/house2.jpg');
%matches=load('/Users/muthuvel/Documents/MATLAB/hw3/data/part2/house_matches.txt');
%Point1=load('/Users/muthuvel/Documents/MATLAB/hw3/data/part2/house1_camera.txt');;
%Point2=load('/Users/muthuvel/Documents/MATLAB/hw3/data/part2/house2_camera.txt');;

a=matches(:,1:2);
b=matches(:,3:4);

%% Calculation of camera centers
% Image 1
[~,~,V]=svd(Point1);
C1=V(:,end);
C1=C1';
% Image 2
[~,~,V]=svd(Point2);
C2=V(:,end);
C2=C2';

%% Triangulation using Linear method
for i=1:1:size(a,1)
    t=a(i,:);
    u=b(i,:);
    P1=Point1;
    P2=Point2;
    A=[(P1(3,:)'*t(2) - P1(2,:)')'; 
       (P1(1,:)'- t(1)*P1(3,:)')';
       (P2(3,:)'*u(2) - P2(2,:)')'; 
       (P2(1,:)'- u(1)*P2(3,:)')'];
    [~,~,V]=svd(A);
    V=V(:,end);
    Val(i,:)=V';
end

end