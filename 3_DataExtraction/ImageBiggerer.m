function [ imageCN ] = ImageBiggerer( image )
%IMAGEBIGGERER Summary of this function goes here
%   image size 108x492 -> size 512x2333

shortSide=512;
longSide=2333; % Used at end just in case ratio is rounded weirdly
%{
ratio1=shortSide/size(image,1);
ratio2=longSide/size(image,2);
ratio = max(ratio1,ratio2);
image=imresize(image, ratio); %%%Upsampling is often quite necessary, you can try it.
if ~isequal(size(image),[shortSide,longSide])
    disp('hai')
    image = image(1:shortSide,1:longSide);
end
size(image)
disp('bai')
%imshow(image);
%pause();

percentile=0.995; % Brightness level for normalization
LImageDim=600; %Real length of the longer FOV side [um]
filtKerSize=75; %Size of the image filtration kernell [um]
filtKerSigm=30; %Sigma of the image filtration kernell [um]

%%%% Precalculations
[RImagePix, CImagePix]=size(image);
LImagePix=max(RImagePix, CImagePix);
pixelation=LImagePix/LImageDim; %Pixels per micron.

image=double(image); %Unfortunately inevitable for image statistics.
filtKer=fspecial('gaussian', round(filtKerSize*pixelation), round(filtKerSigm*pixelation));
background=imfilter(image, filtKer);
normFactor=1./background;
imageC=double(image).*normFactor; %Compensated image

pixels=sort(imageC(:)); %%%% This line and the following one are faster than 
                % "prctile" function from Statistics Toolbox from Mathworks
THR=pixels(round(length(pixels)*percentile)); 
imageCN=uint16(255*imageC/THR); % Compensated and Normalized image
%imshow(imageCN);
%imageCNS=single(imageCN);
%}
imageCN = imresize(image,[shortSide,longSide]) ;

end

