% The function segments the pore space from a cast thin section image.
% Inputs:
    %Image- RGB image of a cast thin section;
    %Hmin- minimum of hue;
    %Hmax- maximum of hue;
    %Smin- minimum of saturation;
    %Vmin- minimum of value.
% Outputs:
    %pores1- the acquired pore space.
function pores1=PoreSegmentation(Image,Hmin,Hmax,Smin,Vmin)
[pp,qq]=size(Image(:,:,1));
ImgHSV = rgb2hsv(Image);% Convert the RGB color model to the HSV color model
[row, col] = ind2sub( size(ImgHSV),...  % segmentation of blue pixels using the HSV color model
    find( ImgHSV(:,:,1)> Hmin & ImgHSV(:,:,1)< Hmax ...
    & ImgHSV(:,:,2)> Smin & ImgHSV(:,:,3)> Vmin ) );
pores=zeros(pp,qq);
for i = 1 : length(row)
    pores(row(i),col(i),:) = 1; % let the value of the pixel in the pore space being 1
end
pores1=logical(pores);
pores1 = bwareaopen(pores1,20,8);  % pore space after deleting the tiny pores
% superposition of segmented pores to the original cast thin section image
Mask_B= 255;  % Image mask
%BP_R, BP_G, and BP_B are the red, green, and blue components of RGB color.
BP_R=Image(:,:,1);BP_G=Image(:,:,2);BP_B=Image(:,:,3);
BP_R(pores1)=Mask_B; BP_G(pores1)=Mask_B; BP_B(pores1)=Mask_B;
Image_H(:,:,1)=BP_R; Image_H(:,:,2)=BP_G; Image_H(:,:,3)=BP_B;
figure
imshow(Image_H);
end