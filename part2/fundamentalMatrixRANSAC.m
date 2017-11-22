%% Function to find matching points and estimate fundamental matrix using RANSAC

function [F,RInlier, matches]=fundamentalMatrixRANSAC(imgLeft,imgRight,flag)
% 10000 iteration for RANASC
%% Inputs
%  imgleft-> image A
%  imgRight-> image B
%  flag-> flag argument to be passed to fundamental_fit()
%  set flag=1 for normalization/ flag=0 for F calculation without normalization

%% Output
%  F->Fundamental Matrix
%  RInlier-> Number of Inliers
%  matches-> 200 putative matches

threshold=2.5; %% earlier used for dist2 calculation, now unused
errorThreshold=5;   %% for RANSAC
r=5;  %%  descriptor size r*r 
r=r+1;
%imgLeft=imread('house1.jpg');
%imgRight=imread('house2.jpg');
newImage=[imgLeft,imgRight];
imgl=imgLeft;
imgr=imgRight;
imgLeft=rgb2gray(imgLeft);
imgRight=rgb2gray(imgRight);
imgLeft=im2double(imgLeft);
imgRight=im2double(imgRight);

%% Feature Extraction from Harris Detector

tic;
fprintf('Extracting corners from Harris detector \n');
[cim_L,row_Leftimage,col_Leftimage]=harris(imgLeft,3,0.001,3,1);  %%Thresh = 0.02 ideal as it detects lot of features
[cim_R,row_Rightimage,col_Rightimage]=harris(imgRight,3,0.001,3,1);
toc;

close all;
tic
fprintf('Extracting decriptors for 1st Image \n');
imgLeft=padarray(imgLeft,[r,r],'replicate','both');
imgRight=padarray(imgRight,[r,r],'replicate','both');

%% Extracting descriptors from corners

wait=waitbar(0,'Extracting decriptors for 2nd Image');
for i=r:1:length(col_Leftimage)-r
    idxi=row_Leftimage(i)+r;
    idxj=col_Leftimage(i)+r;
    
    descriptorLeft{i}=imgLeft(idxi-r:idxi+r,idxj-r:idxj+r);
    
    descriptorLeft{i}=(descriptorLeft{i}(:))';
    
    if i==r
        descpl=descriptorLeft{i};
    end
    if i~=r
        descpl=cat(1,descpl,descriptorLeft{i});
    end
    waitbar(i/length(col_Leftimage),wait);
    
    
end
delete(wait);



wait=waitbar(0,'Extracting decriptors');
for i=r:1:length(col_Rightimage)-r
    idxi=row_Rightimage(i)+r;
    idxj=col_Rightimage(i)+r;
    
    descriptorRight{i}=imgRight(idxi-r:idxi+r,idxj-r:idxj+r);
    
    descriptorRight{i}=(descriptorRight{i}(:))';
    if i==r
        descpr=descriptorRight{i};
    end
    if i~=r
        descpr=cat(1,descpr,descriptorRight{i});
    end
    waitbar(i/length(col_Rightimage),wait);
end
clear i;
clear j;
toc;
delete(wait);

%% Pairwise difference calculation

fprintf('Pairwise distance calculation \n');
tic;
%wait=waitbar(0,'Pairwise distance calculation');
pairwiseDiff=zeros(1,10);

for i=r:1:length(descriptorLeft)-r
    for j=r:1:length(descriptorRight)-r
        pairwiseDiff(i,j)=dist2(descriptorLeft{i},descriptorRight{j});
    end
end

pairwiseDiff(:,1:r-1)=1000000; 
pairwiseDiff(1:r-1,:)=1000000;%% some big value
pairwiseDiff(length(descriptorLeft)-r:length(descriptorLeft),:)=10000000;
pairwiseDiff(:,length(descriptorRight)-r:length(descriptorRight))=10000000;
toc;
clear descriptorLeft;
clear descriptorRight;
clear i;
clear j
clear idxi;
clear idxj;
rowL=[];
colL=[];
rowR=[];
colR=[];
v=[];
wait=waitbar(0,'Extracting putative matches');

for i=1:1:200
    a=min(min(pairwiseDiff));
    [idxi,idxj]=find(pairwiseDiff==a);
    rowL=[rowL; row_Leftimage(idxi)];
    colL=[colL; col_Leftimage(idxi)];
    rowR=[rowR; row_Rightimage(idxj)];
    colR=[colR; col_Rightimage(idxj)];
    v=[v,a];
    pairwiseDiff(idxi,:)=100000; %%some big value;
    pairwiseDiff(:,idxj)=100000;
    
    waitbar(i/200,wait)
end

delete(wait);
clear idxi;
clear idxj;
clear row_pdl;
clear row_pdr;
clear col_pdl;
clear col_pdr;

%% RANSAC for estimation of fundamental matrix
x=rowL;
y=colL;
X=rowR;
Y=colR;
threshold=errorThreshold;
fprintf('Calculating Inliers and Ouliers\n')
tic;
wait=waitbar(0,'Calculating Inliers and Ouliers');

for iteration=1:1:10000
    x1=[];
    X1=[];
    y1=[];
    Y1=[];
    idx=randsample(length(X),8);
    for j=1:1:8
        x1=[x1;x(idx(j))];
        y1=[y1;y(idx(j))];
        X1=[X1;X(idx(j))];
        Y1=[Y1;Y(idx(j))];
    end
    match_1=[x1,y1,X1,Y1];
    
    % Calls fundamental_fit function for finding F from matching points
    % detected
    F=fundamental_fit(match_1,flag);
    
    Inlier=0;
    Outlier=0;
    %for j=1:1:length(X)
        %error(j)=[X(j), Y(j), 1]*F*[x(j);y(j);1];
        
       % if (abs(error(j))<threshold)
       %     Inlier=Inlier+1;
       %     continue;
       % end
       % Outlier=Outlier+1;
    %end

    L = (F * [[x,y] ones(length(x),1)]')'; 
    L = L ./ repmat(sqrt(L(:,1).^2 + L(:,2).^2), 1, 3); 
    pt_line_dist = sum(L .* [[X,Y] ones(length(X),1)],2);
    %closest_pt = matches(:,3:4) - L(:,1:2) .* repmat(pt_line_dist, 1, 2);
    
    Inlier=size(find(abs(pt_line_dist)<errorThreshold));
    coor=find(abs(pt_line_dist)<errorThreshold);
    Outlier=size(find(abs(pt_line_dist)>=errorThreshold));
    
    if (iteration == 1)
        RInlier=Inlier(1);
        ROutlier=Outlier(1);
        RF=F;
        RCOOR=coor;
    end
    
    if (Inlier(1)>RInlier)
        RInlier=Inlier(1);
        ROutlier=Outlier(1);
        RF=F;
        RCOOR=coor;
    end
    waitbar(iteration/10000,wait);
        
    
end

F=RF; %% Assigning output value
matches=[rowL,colL,rowR,colR];
matches=matches(RCOOR,:);
toc;
delete(wait);
clear idx;
clear Inlier;
clear Outlier;
clear H;
clear X;
clear X;
clear Y;
clear P;
clear h;
clear Z;
clear mat_y;
clear mat_x;
clear iteration;
clear InlierLeft;
clear InlierRIght;
clear IR;
clear IL;
clear v;
clear V;
clear Vector;
clear wait;
clear i;
clear idxi;
clear A;
clear B;


end




