%% Function to estimate homography
%% Uses RANSAC Algorithm - 1000 iterations
%% Change input image in line 14 and 15
%% De-comment line 141 and 142 to visualize putative matches 
clc;
clear;
close all;

%threshold=2.5; %% for dist 2 NOT USED NOW
SSD_Threshold=50;   %% for RANSAC
r=5;  %% for descriptor size r*r 
r=r+1;

imgLeft=imread ('uttower/left.jpg');
imgRight=imread ('uttower/right.jpg');
newImage=[imgLeft,imgRight];
imgl=imgLeft;
imgr=imgRight;
imgLeft=rgb2gray(imgLeft);
imgRight=rgb2gray(imgRight);
imgLeft=im2double(imgLeft);
imgRight=im2double(imgRight);

%% Extraction of corners from Harris Corner

tic
fprintf('Extracting corners from Harris detector \n');
[cim_L,row_Leftimage,col_Leftimage]=harris(imgLeft,3,0.05,3,1);  %%Thresh = 0.02 ideal as it detects lot of features
[cim_R,row_Rightimage,col_Rightimage]=harris(imgRight,3,0.05,3,1);
toc;

%close all;
tic
fprintf('Extracting decriptors \n');
imgLeft=padarray(imgLeft,[r,r],'replicate','both');
imgRight=padarray(imgRight,[r,r],'replicate','both');
%match_plot(imgl,imgr,[col_Leftimage,row_Leftimage],[col_Rightimage,row_Rightimage]);

%% Extracting Descriptors from corners

wait=waitbar(0,'Extracting decriptors');
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

%% Pairwise distance calculation from dist2

fprintf('Pairwise distance calculation \n');
tic;
% wait=waitbar(0,'Pairwise distance calculation');
pairwiseDiff=zeros(1,10);

for i=r:1:length(descriptorLeft)-r
    for j=r:1:length(descriptorRight)-r
        pairwiseDiff(i,j)=dist2(descriptorLeft{i},descriptorRight{j});
    end
end

pairwiseDiff(:,1:r-1)=1000000; 
pairwiseDiff(1:r-1,:)=1000000;%% some big valuees to eliminate recurrence of same points
pairwiseDiff(length(descriptorLeft)-r:length(descriptorLeft),:)=10000000;
pairwiseDiff(:,length(descriptorRight)-r:length(descriptorRight))=10000000;

% delete(wait);
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

% rowL=row_pdl(2:length(row_pdl));  NOT USED
% colL=col_pdl(2:length(col_pdl));
% rowR=row_pdr(2:length(row_pdr));
% colR=col_pdr(2:length(col_pdr));
% v=v(2:length(v));
% Pointleft=[colL;rowL]';  USED FOR VERIFICATION
% PointRight=[colR;rowR]';

%ts=rowR;
%rowR=colR;
%colR=ts;

%ts=colR;
%colL=rowL;
%rowL=ts;

%match_plot(imgl,imgr,[colL,rowL],[colR, rowR]);
%title('Putative matches');


clear idxi;
clear idxj;
clear row_pdl;
clear row_pdr;
clear col_pdl;
clear col_pdr;

%% RANSAC for homography estimation

fprintf('Calculating Inliers and Ouliers\n')
tic;
wait=waitbar(0,'Calculating Inliers and Ouliers');

for iteration=1:1:10000
    idx=randsample(length(v),4);
    P=[];
    for i=1:1:4
        idxi=idx(i);
        x=rowL(idxi);
        y=colL(idxi);
        X=rowR(idxi);
        Y=colR(idxi);
   
        P=[P; -x, -y, -1, 0, 0, 0, x*X, y*X, X];
        P=[P; 0, 0, 0, -x, -y, -1, x*Y, y*Y, Y];  
    end
    % Calculation of homography matrix
    [A,B,V]=svd(P);
    h=V(:,end);
    h=h/h(9);
    H=reshape(h,[3,3]);
    H=H';
    Hom{iteration}=H;
    
    X=[];
    Y=[];
    Z=[];  
    
    for i=1:1:length(colR)
        Points=[rowL(i); colL(i); 1];
        [Vector]=H*Points;
        Points=Points./0.0001;
        X=[X,Vector(1)];
        Y=[Y,Vector(2)];
        Z=[Z,Vector(3)];
    end
    
    clear i;

    Inlier=0;
    Outlier=0;
    Irowl=[];
    Icoll=[];
    Irowr=[];
    Icolr=[]; 
    InlierRight=[];
    
    for i=1:1:length(colL) %%%%% try chagnging rowR bloe to colR and likewise
       SSD(i)=sqrt(((X(i)-rowR(i))^2)+((Y(i)-colR(i))^2));
       if (SSD(i)<=SSD_Threshold)
           Inlier=1+Inlier;
           Irowl=[Irowl,rowL(i)];
           Icoll=[Icoll,colL(i)];
           Irowr=[Irowr,rowR(i)];
           Icolr=[Icolr,colR(i)];
       end
       if (SSD(i)>SSD_Threshold)
           Outlier=1+Outlier;
       end
    end
    Inl(iteration)=Inlier;
    if(iteration==1)
            RANSAC_INLIER=Inlier;
            RANSAC_OUTLIER=Outlier;
            RANSAC_Coordinates_for_homography=idx;
            RH=H;
            RIrowl=Irowl;
            RIcoll=Icoll;
            RIrowr=Irowr;
            RIcolr=Icolr;
    end
    if (iteration~=1)
        if (RANSAC_INLIER<Inlier)
            RANSAC_INLIER=Inlier;
            RANSAC_OUTLIER=Outlier;
            RANSAC_Coordinates_for_homography=idx;
            RH=H;
            RIrowl=Irowl;
            RIcoll=Icoll;
            RIrowr=Irowr;
            RIcolr=Icolr;
        end
    end
    waitbar(iteration/10000,wait);   
end
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
clear Icoll;
clear Icolr;
clear Irowl;
clear Irowr;
clear imgRight;
clear imgLeft;
clear i;
clear B;
clear t;
clear Points;

fprintf("Number of Inliers:");disp(RANSAC_INLIER);
fprintf("Residue=");disp((sum(SSD)/200));

%% Remove this
Pointleft=[RIcoll;RIrowl;]';
PointRight=[RIcolr;RIrowr;]';
RH=RH';

match_plot(imgl,imgr,Pointleft,PointRight);
title('Lines joining inliers');

figure(4)
imshow(imgl);
hold on;
plot(Pointleft(:,1),Pointleft(:,2),'sr');
title('Inlier points in image 1');
figure(5)
imshow(imgr);
hold on;
plot(PointRight(:,1),PointRight(:,2),'sr');
title('Inlier points in image 2');

%% Make transform and image stitching


T=maketform('projective',RH);
[B,xdata,ydata] = imtransform(imgr,T);
figure(6);
imshow(B);
title('Warped image');

% Code for Stitching
%X=[min(1,xdata(1)) max(size(imgl,2), xdata(2))];
%Y=[min(1,ydata(1)) max(size(imgl,1), ydata(2))];
%img1 = imtransform(imgr,T);
%img2 = imtransform(imgl, maketform('affine',eye(3)));
%finalImage = img1 + img2;
%imshow(finalImage)
