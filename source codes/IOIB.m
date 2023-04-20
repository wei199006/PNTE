% The function figures out inlets, outlets, and intersectional points & blocks of the pore network.
% Inputs:
    % pp- the row number of the cast thin section image;
    % qq- the column number of the cast thin section image;
    % sp1- the pore network.
% Outputs:
    % bw- the final pore network image;
    % pois- the collections of inlets, outlets, and intersectional points & blocks of the pore network;
    % sp- the inlets of the pore network;
    % ep- the outlets of the pore network;
    % mp- the intersectional points of the pore network;
    % bp- the intersectional blocks of the pore network.
function [bw,pois,sp,ep,mp,bp]=IOIB(pp,qq,sp1)
% delete the noise pixels mistaken identified as inputs or outputs
for i= 1:pp
    lsp= [0;sp1(:,2);0];
    rsp= [0;sp1(:,qq-1);0];
    if sp1(i,1)==1 && lsp(i)+lsp(i+1)+lsp(i+2)<0.5
       sp1(i,1)= 0; 
    end
    if sp1(i,qq)==1 && rsp(i)+rsp(i+1)+rsp(i+2)<0.5
       sp1(i,qq)= 0; 
    end
end
sp1= bwareaopen(sp1,50,8);% delete short pore networks that are not connected to the main pore network.
bw= sp1;
[mm,nn]= size(bw);
% Identification for inlets
ii= 1;
for i= 1:mm
    if (bw(i,1)== 1)
        sp(ii,1)= i; % position vector for inlets
        sp(ii,2)= 1;
        ii=ii+1;
    end
end
% Identification for outlets
ii= 1;
for i= 1:mm
    if (bw(i,nn)== 1)
        ep(ii,1) = i; % position vector for outlets
        ep(ii,2) = nn;
        ii= ii+1;
    end
end
% Identification for pore network intersectional points
ii=1;
for i=2:mm-1
    for j=2:nn-1
        la=[bw(i-1,j), bw(i-1,j+1), bw(i,j+1), bw(i+1,j+1),...
            bw(i+1,j), bw(i+1,j-1), bw(i,j-1), bw(i-1,j-1)];
        sumv= sum(la);
        for k=2:7
            if(la(k)>0.5 && la(k-1)+la(k)+la(k+1)>1.5)
                sumv=0;
            end
        end
        if(la(8)>0.5 && la(7)+la(8)+la(1)>1.5)
            sumv=0;
        end
        if(la(1)>0.5 && la(8)+la(1)+la(2)>1.5)
            sumv=0;
        end
        if (bw(i,j)>0.5 && sumv>2.5) % the value of it equals 1 & at least, three of its neiber points equal 1
            mp(ii,1)= i; % position vector for intersections
            mp(ii,2)= j;
            ii=ii+1;
        end
    end
end
% Identification for pore network intersectional blocks
blocksn= 1;%number of  pixels identified as points of intersectional blocks
for i=1:pp-1
    for j=1:qq-1
        if abs(bw(i,j)+bw(i,j+1)+bw(i+1,j)+bw(i+1,j+1)-4)<10^(-6)
            blocksxy(blocksn,1)= i;%x coordinate of  pixels identified as points of intersectional blocks
            blocksxy(blocksn,2)= j;%y coordinate of  pixels identified as points of intersectional blocks
            blocksn= blocksn+1;
        end
    end
end
blocksn= 1;
if exist('blocksxy','var')
    for i=1:length(blocksxy)
        blocks1{blocksn}= [blocksxy(i,1) blocksxy(i,2);blocksxy(i,1) blocksxy(i,2)+1;...
            blocksxy(i,1)+1 blocksxy(i,2);blocksxy(i,1)+1 blocksxy(i,2)+1 ];
        blocksn= blocksn+1;
    end
    blocksn= 1;
    blocks{1}=blocks1{1};
    for i=2:length(blocks1)
        lj=length(blocks);
        flag=1;
        for j=1:lj
            % if the intersection of two blocks is not empty then unit the two blocks
            if ~isempty(intersect(blocks{j},blocks1{i},'rows'))
                blocks{blocksn}=union(blocks{j},blocks1{i},'rows');
                flag= 0;
                break;
            end
        end
        if abs(flag-1)<10^(-6)
            blocksn=blocksn+1;
            blocks{blocksn}=blocks1{i};
        end
    end
    for i=1:length(blocks)
        bp(i,1)= floor(mean(blocks{i}(:,1))); % position vector for intersectional blocks
        bp(i,2)= floor(mean(blocks{i}(:,2))); 
    end
    pois=[sp; ep; mp; bp]; % pois- a vector containing all the vertex
else
    pois=[sp; ep; mp]; % pois- a vector containing all the vertex
end
% figure plot
% figure
% imshow(~bw);
% hold on 
% scatter(sp(:,2),sp(:,1),'green','square','filled'); % plot inlets
% scatter(ep(:,2),ep(:,1),'green','diamond','filled'); % plot outlets
% scatter(mp(:,2),mp(:,1),'green','o','filled'); % plot intersectional points
% if exist('blocksxy','var')
%     scatter(bp(:,2),bp(:,1),'blue','o','filled'); % plot intersectional blocks
% end
% [tm,tn]=size(pois);
% % Mark the numbers on the intersections
% for i=1:tm
%     c = num2str(i);
%     text(pois(i,2)+5,pois(i,1),c);
% end
% hold off
end