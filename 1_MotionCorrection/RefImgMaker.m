function [ reference ] = RefImgMaker( stack )
% Finds average of a stack of images

imageN=size(stack,3);
imageN = 75
% Should be between 0.5 and 1.5, the lower it is the brighter the image
meanFactor = imageN*0.3; 
% meanFactor=15
meanStack = stack(:,:,1)/meanFactor;
for i=2:imageN
% for i=190:250
    meanStack = meanStack + stack(:,:,i)/meanFactor;
end
%  reference = uint16(meanStack);
reference = imresize(meanStack,10);

imshow(reference,'Border','tight');
saveas(gca,'DG Reference.png');

end