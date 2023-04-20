% Copyright Wei Li 2023/01/22.
% Please report any bug to liwei2014@email.cugb.edu.cn or weili.rockphysics@gmail.com.
% Referred manuscript title: "A program for estimating the tortuosity of rocks from cast thin-section images".
% The PNTE algorithm is short for the pore network tortuosity estimation.
% The algorithm provides an estimation of the tortuosity for a pore network of a cast thin section image.
clc;
clear;
close all;
Image=imread('test image.jpg');
%% Input parameters
% Tortuosity definition flags
dtype= 1; %=1 width of the image, =2 straight-line distance between the inlet and the outlet.
ttype= 2; %=1 average tortuosity among all pi, =2 minimum tortuosity among all pi.
% Pore space segmentation parameters (HSV color model)
Hmin= 0.6; % minimum of hue.
Hmax= 0.7; % maximum of hue.
Smin= 0.3; % minimum of saturation.
Vmin= 0.2; % minimum of value.
% Dilation factor Df
dilaf= 35; % dilaf- dilation factor Df (should be a multiple of 5).
details= 0; %=1 show the figures referred to the gradual revisions of the pore network, =0 skip the referred figures.
%% Pore space segmentation (Procedure A)
[pp,qq]= size(Image(:,:,1));% pp and qq - the row and the column numbers of the cast thin section image, respectively.
pores1= PoreSegmentation(Image,Hmin,Hmax,Smin,Vmin);
%% Pore network generation (Procedure B)
sp1= zeros(pp,qq); % sp1-the pore network.
sp1= pathrevision(pp,qq,dilaf,pores1,sp1,details); % the pore network acquired before the revision.
%% Pore network revision (Procedure C1-C5)
dilafv= dilaf:-5:5; % gradually reduce the dilation factor.
for i= 1:length(dilafv)
    sp1= pathrevision(pp,qq,dilafv(i),pores1,sp1,details); % the pore network acquired after the revision.
end
%% Removement of complex branches (Procedure C6)
dilaf= 3;
zpore= zeros(pp,qq);%empty matrix
sp1= pathrevision(pp,qq,dilaf,sp1,zpore,details); % the pore network acquired before the revision
%% Identifying the inlets, outlets, and intersections of the pore network (Procedure D1 & D2)
[bw,pois,sp,ep,mp,bp]= IOIB(pp,qq,sp1);
[mm,nn]= size(bw);
%% generate distance-weighted graph (Procedure D3)
[dis,tnote,tpath,tbns]= DWG(bw,pois);
%% Floyd-Warshall Algorithm (Procedure E)
[dis,tnote]= FWalgorithm(dis,tnote,pois);
%% Tortuosity estimation (Procedure F)
[taupi,taun,pi,nt]= tau(sp,ep,dis,dtype,ttype,qq);
%% Trace the shortest paths
[nettrace,netpath]= SPtrace(mm,sp,pi,bw,tnote,tpath,tbns,nt);
%% Figure Plot
% Pore network image marked with all the vertices
figure
imshow(~bw);
hold on
bpmp= setdiff(pois,union(sp,ep,'rows'),'rows');
scatter(bpmp(:,2),bpmp(:,1),'green','o','filled');
scatter(sp(:,2),sp(:,1),'red','square','filled');
scatter(ep(:,2),ep(:,1),'blue','diamond','filled');
[tm,tn]=size(pois);
for i=1:tm
    c= num2str(i); % convert the numbers to chars
    text(pois(i,2)+5,pois(i,1),c); % mark the numbers of the vertices on the figure
end
% plot all pi
for i= 1:length(netpath)
    if ~isnan(netpath{i})
        pathplot=netpath{i};
        plot(pathplot(:,2),pathplot(:,1),'g');
    end
end
hold off

% Superposition of the pore network and all pi onto the original cast thin section image
figure
% plot the pore network on the original cast thin-section image
bwplot= bw;
Mask_B= 255; % image mask
BP_R= Image(:,:,1);BP_G= Image(:,:,2);BP_B= Image(:,:,3);
BP_R(bwplot)= Mask_B; BP_G(bwplot)= Mask_B; BP_B(bwplot)= Mask_B;
Image_L(:,:,1)= BP_R; Image_L(:,:,2)= BP_G; Image_L(:,:,3)= BP_B;
imshow(Image_L);
hold on
% plot all pi
for i= 1:length(netpath)
    if ~isnan(netpath{i})
        pathplot= netpath{i};
        plot(pathplot(:,2),pathplot(:,1),'g','Linewidth', 1);
    end
end
% plot the shortest path among all pi
pathplot= netpath{nt};
plot(pathplot(:,2),pathplot(:,1),'r','Linewidth', 1);
hold off