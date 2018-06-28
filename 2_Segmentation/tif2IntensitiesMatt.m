function tif2IntensitiesMatt(stack, maxDisplacement,outputName, region)

%region = 'DG s'; % Which region of the brain is this movie in? (DG, CA1)
               % used in segmCorr to determine cell sizes
               % Pixels per micron: 3.8883 (/2 ?)
shortSide=512;
%percentileBrightness=95; %How to set the threshold in reference seeking process.
neuropilCorr=0.2; %How many times should be neuropil subtracted. Chen2013 used 0.3;
neuropilDistance=1; %How many times should be the sigma of considered neuropil area
neuropilNPixel=1000; %%%From how many randomly distributed neuropil point the mean neuropil signal should be measured.
dilationParameter=5; %%%Roughly... how far (in pixels) should be the enabled area for picking the neuronpil pixels.

fileName= outputName;%input('Type the name for your data file.', 's');
%fileName='TSeriesIntensities001.mat'; %%%The food for TwoPhotonProcessor (our software published in 2013)
%%%%I'will send you a new version of TPP on request. 

% stack(1:maxDisplacement,:,:)=[];
% stack(:,1:maxDisplacement,:)=[];
% stack(end-maxDisplacement+1:end,:,:)=[];
% stack(:,end-maxDisplacement+1:end,:)=[];
[maskStructInt, imageCN, hText, rgbBlueCO, maskStack]=tiffIntegratorMaskM_Prague(stack,shortSide,neuropilCorr, neuropilDistance,neuropilNPixel,dilationParameter,fileName,region);


neuronN=size(maskStructInt,2);
numOfImages=size(maskStructInt(1,1).Intensity,2);
intavg=zeros(neuronN,numOfImages);
intavgNP=zeros(neuronN,numOfImages);
intavgNPCorr=zeros(neuronN,numOfImages);
%%%This is the actual calcium-related fluorescence trace.
positionTimeCorr=zeros(neuronN,1); %%%This variable is useful to refine the time resolution. Even if you sample with a certain
%%%fps, you still know that the image is scanned from top to bottom-> Then
%%%the precise time when a neuron is recorded is
%%%NthFrame/fps+(yCoorNeuron/SizeOfFrame)/fps. It may improve the precision at least 10 times and may prove useful.
for index=1:neuronN    
    intavg(index, :)=maskStructInt(1,index).Intensity;
    intavgNP(index, :)=maskStructInt(1,index).IntensityNP;
    intavgNPCorr(index, :)=maskStructInt(1,index).IntensityNPCorr;
    positionTimeCorr(index)=maskStructInt(1,index).CenterCoor(1)/shortSide;
end

%rgb = zeros(size(rgbBlueCO,1),size(rgbBlueCO,2),size(stack,3),3);
for index=1:size(stack,3)
    stack(1,1,index);
    stack(:,:,index) = stack(:,:,index) *3;%rgbBlueCO;
end
%toTiff(stack, [fileName,' CUT outlined']);

fileName(end+1:end+4)='.mat';
save(fileName,'intavg','positionTimeCorr', 'intavgNP', 'intavgNPCorr', 'maskStructInt', 'imageCN', 'numOfImages', 'hText', 'maskStack');

end