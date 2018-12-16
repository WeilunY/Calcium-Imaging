function [maskStructInt, imageCN, hText, rgbBlueCO, maskStack]=tiffIntegratorMaskM_Prague(stack,shortSide,neuropilCorr, neuropilDistance, neuropilNPixel,dilationParameter, fileName, region)

meanTiff=mean(stack(:,:,1:200),3);%meanTiff=mean(stack(:,:,1:200),3); %Seeking reference for segmentation
ratio=shortSide/size(meanTiff,1);
meanTiff=imresize(meanTiff, ratio); %%%Upsampling is often quite necessary, you can try it.

[maskStructInt,imageCN, outerSpace, neuronR, rgbBlueCO, maskStack]=segmCorr(meanTiff, stack, region, fileName); %%%The core
points=find(outerSpace<=0);
outerNew=zeros(size(outerSpace));
outerNew(points)=1;
a=ones(dilationParameter); %%%%%tady tuto 5ku navázat na rozlišení 
outerSpace=imdilate(outerNew,a);
imageSize=size(meanTiff);

% imshow(imageCN);figure(gcf);
% pause;

numImages=size(stack,3);
maskStructInt(1,1).Intensity=zeros(1,numImages);
maskStructInt(1,1).IntensityNP=zeros(1,numImages);
maskStructInt(1,1).IntensityNPCorr=zeros(1,numImages);
maskStructInt(1,1).IntensityWC=zeros(1,numImages);

neuronN=size(maskStructInt,2);
maskStructInt(1,1).NPixList=zeros(neuropilNPixel,1); %%%Indices of corresponding neuropil pixels
h = waitbar(0,'Picking the proximal neuropil pixels');
for indexN=1:neuronN %%%%%Zde se pro každý neuron natrvalo napoèítají neuropilové pixely (v indexovane reprezentaci v celem obrazku)
    waitbar(indexN / neuronN) 
    pixelNumber=1;
    middlePoint=round(maskStructInt(1,indexN).CenterCoor); %%%Najdu centrum neuronu, coor(1) je øádek, coor(2) je sloupec
    maskStructInt(1,indexN).NPixList(neuropilNPixel)=0;
while maskStructInt(1,indexN).NPixList(neuropilNPixel)==0
    %%%zkoušim body v polárních souøadnicích, polomìr je randn, smìr isotropní
        pixDist=neuropilDistance*neuronR*abs(randn);
        pixAngle=2*pi*rand;
        [X,Y] = pol2cart(pixAngle,pixDist);
        newPoint=[middlePoint(1)+round(Y) middlePoint(2)+round(X)];
    if newPoint(1)>=1 && newPoint(1)<=imageSize(1) && newPoint(2)>=1 && newPoint(2)<=imageSize(2) && outerSpace(newPoint(1),newPoint(2))==0  
        maskStructInt(1,indexN).NPixList(pixelNumber)= sub2ind(size(meanTiff),newPoint(1),newPoint(2));
        pixelNumber=pixelNumber+1;
    end
end
% imshow(outerSpace);
% hold on
% [I,J] = ind2sub(size(meanTiff),maskStructInt(1,indexN).NPixList);
% plot(J,I, '.');
% pause;
% close;
end

close(h);

for indexMarker=1:neuronN %%%Here it produces a nice annotated reference image
    hText = text(round(maskStructInt(1,indexMarker).CenterCoor(2)),round(maskStructInt(1,indexMarker).CenterCoor(1)),num2str(indexMarker),'Color',[1 0 0],'FontSize',6);
end
%pause();
%hText   %%% This has all of the cell labelling
%pause();

h = waitbar(0,'Calculating the intensities...');
for k = 1:numImages %%%For cycle to go through individual images
    waitbar(k / numImages)        
    imgTMP = stack(:,:,k);
    imgTMP=imresize(imgTMP, ratio);
    for indexN=1:neuronN %%%For cycle to obtain average intensities within each image
            %%%
        middlePoint=round(maskStructInt(1,indexN).CenterCoor);
        halfNeuron=(maskStructInt(1,indexN).surroundSize(1)-1)/2;
        dataN=imgTMP(middlePoint(1)-halfNeuron:middlePoint(1)+halfNeuron, middlePoint(2)-halfNeuron:middlePoint(2)+halfNeuron);

        maskStructInt(1,indexN).Intensity(k)=mean(dataN(maskStructInt(1,indexN).ringPixList)); %%% Intensity from the neuron body
        maskStructInt(1,indexN).IntensityNP(k)=mean(imgTMP(maskStructInt(1,indexN).NPixList)); %%% Intensity from the surrounding neuropil
        maskStructInt(1,indexN).IntensityWC(k)=mean(dataN(maskStructInt(1,indexN).WCPixList)); 
                                                                                               %%% Body intensity corrected for neuropil.
        maskStructInt(1,indexN).IntensityNPCorr(k)=maskStructInt(1,indexN).Intensity(k)-neuropilCorr*maskStructInt(1,indexN).IntensityNP(k); 
    end
end
close(h);

fileName(end+1:end+4)='.fig';
saveas(gca,fileName)
fileName(end-3:end)='.png';
set(gcf,'PaperPositionMode','auto')
print(gcf, '-r300', '-dpng',fileName);

end
