% The function is designed to generate and revise the pore network from the pore space.
% Inputs:
    %pp- the number of the rows of the cast thin section image;
    %qq- the number of the coloumns of the cast thin section image;
    %dilaf- dilation factor;
    %pores1- the pore space;
    %sp1- the pore network;
    %details- flag parameter defined in the main function.
% Outputs:
    %osp1- the new pore network.
% Notes: The procedures referred to the generation and the revision of the pore network are similar.
% For consicely, both procedures are achieved through this function.
% In the generation, the input spl is an empty matrix (see procedure B).
% In the revision, the input spl is the pore network that needs revision (see procedure C1-C5).
function osp1=pathrevision(pp,qq,dilaf,pores1,sp1,details)
% image padding (Procedure B1)
nfig= zeros(pp+dilaf*4,qq+dilaf*4); %nfig- pore space image after padding.
nfig(dilaf*2+1:dilaf*2+pp,dilaf*2+1:dilaf*2+qq)= pores1;
nfig1= zeros(pp+dilaf*4,qq+dilaf*4); %nfig1- pore network image after padding.
nfig1(dilaf*2+1:dilaf*2+pp,dilaf*2+1:dilaf*2+qq)= sp1;
se = strel('disk',dilaf); %Circular operator
dnfig = imdilate(nfig,se); %dilation (Procedure B2)
[emm,enn]=size(nfig);
udnfig=zeros(emm,enn);
% unit the pore space and the pore network
for i=1:pp+dilaf*4
    for j=1:qq+dilaf*4
        if (nfig1(i,j)==1 || dnfig(i,j)==1)
            udnfig(i,j)= 1; % udnfig-union of pore space and pore network
        end
    end
end
if abs(details-1)<10^(-10)
    figure;
    imshow(~udnfig);
    hold on
    %add a black border to the padded image
    plot([1,enn],[1,1],'b');
    plot([1,1],[1,emm],'b');
    plot([enn,enn],[1,emm],'b');
    plot([1,enn],[emm,emm],'b');
    %add a red border to mark the original size of the image (before padding)
    plot([dilaf*2+1,dilaf*2+qq],[dilaf*2+1,dilaf*2+1],'r');
    plot([dilaf*2+1,dilaf*2+1],[dilaf*2+1,dilaf*2+pp],'r');
    plot([dilaf*2+qq,dilaf*2+qq],[dilaf*2+1,dilaf*2+pp],'r');
    plot([dilaf*2+1,dilaf*2+qq],[dilaf*2+pp,dilaf*2+pp],'r');
    hold off
end
%skeletonize the pore space (Procedure B3)
nsk = bwmorph(udnfig,'skel',Inf); %nsk- the pore skeleton image
%remove image padding (Procedure B4)
nskel=nsk(dilaf*2+1:dilaf*2+pp,dilaf*2+1:dilaf*2+qq); %nskel- the pore skeleton image without padding
%prune the branches (Procedure B5)
for i = 1:100
    nprun = bwmorph(nskel,'spur');
    nskel = nprun;
end
osp1= nskel; %ospl- the pore network image
end