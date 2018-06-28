function [maskStructure,imageCN, outerSpace, neuronR, rgbBlueCO, maskStack]=segmCorr(image, stack, region, fileName)
% region = DG or CA1 or...

%%% Parameters
percentile=0.995; % Brightness level for normalization
LImageDim=600; % Real length of the longer FOV side [um]
filtKerSize=75; % Size of the image filtration kernell [um]
filtKerSigm=30; % Sigma of the image filtration kernell [um]
neuronSizeOut=12.5; % Body diameter of a typical neuron representing given brain area; e.g.in L2/3 cortex it is 12-15um,
neuronSizeIn=8.5; % Nucleus diameter of a typical neuron representing given brain area; e.g. in L2/3 cortex it is 8-10um
mode=2; % mode=1 for forced circular neurons (faster, robust), mode=2 for irregular neurons with some constraints 
                                                        %(slower, possible bugs, but theoretically most accurate)
%inflateParam=30;   % v hodnotì procent od støedu neuronu ke kraji
%blurDiam=5; %%% Parameter of disc blurring
%inflateCoeff=inflateParam/100+1;%%% ????
scoreTHR=0.15; %%%Threshold for the inner findDonuts scoring function. Please modify. This needs to be adjusted according to your data. 
                 %The higher the scoreTHR, the lower sensitivity, but the higher specifity. Vice versa.
springCoeff=0.3; %%%After a putative neuron is transformed into polar coordinates. There is a routine to find its ratio(~angle) curve.
%%%This coefficient sets the stiffness against deformation from pure circle (the deformation is pushed by the maximal intensity seeking). 
qualityCheckHandPicked=1; %%% Equals 1 if you want hand picked cells to undergo quality checking
                          %%% Equals 0 if no quality check
              %%% Very Small cells removed by default in findDonuts.m line 125
% Change neuron sizes (in microns), Pixels per micron: 3.8883
%%% Readjusting parameters for CA1
if strcmp(region,'CA1')
    neuronSizeOut = neuronSizeOut+2.5;
    neuronSizeIn = neuronSizeIn+2.5;
    filtKerSize = 75; %Size of the image filtration kernell [um]
    filtKerSigm = 30; %Sigma of the image filtration kernell [um]
elseif strcmp(region,'CA3')
    neuronSizeOut = neuronSizeOut+5.5;
    neuronSizeIn = neuronSizeIn+5.5;
elseif strcmp(region,'Hilus')
    
elseif strcmp(region,'DG')
    neuronSizeOut = neuronSizeOut+1;
    neuronSizeIn = neuronSizeIn+1;

elseif strcmp(region,'DG l')
    neuronSizeOut = neuronSizeOut+1.5;
    neuronSizeIn = neuronSizeIn+1.5;
     
elseif strcmp(region,'DG s') % DG but larger field of view
    neuronSizeOut = neuronSizeOut;
    neuronSizeIn = neuronSizeIn;
    
elseif strcmp(region,'Test')
     neuronSizeOut = 14;
     neuronSizeIn = 12;
end
%neuronSizeOut = 30
%neuronSizeIn = 1
              
%%%% Precalculations
[RImagePix, CImagePix]=size(image);
LImagePix=max(RImagePix, CImagePix);
pixelation=LImagePix/LImageDim; %Pixels per micron.

[rr cc] = meshgrid(1:round(neuronSizeOut*pixelation)); %% Proper GCaMP pattern
kernell=and(round(neuronSizeIn*pixelation/2)<= sqrt((rr-round(neuronSizeOut*pixelation/2)).^2+...
    (cc-round(neuronSizeOut*pixelation/2)).^2), sqrt((rr-round(neuronSizeOut*pixelation/2)).^2+...
    (cc-round(neuronSizeOut*pixelation/2)).^2)<=round(neuronSizeOut*pixelation/2));
if isempty(find(kernell(:,end))) %%%Cropping out the empty column/row (due to rounding) 
kernell(:,end)=[];
kernell(end,:)=[];
end
%imagesc(kernell);
%pause;
%%%

%%% Image preparation
image=double(image); %Unfortunately inevitable for image statistics.
filtKer=fspecial('gaussian', round(filtKerSize*pixelation), round(filtKerSigm*pixelation));
background=imfilter(image, filtKer);
normFactor=1./background;
imageC=double(image).*normFactor; %Compensated image

pixels=sort(imageC(:)); %%%% This line and the following one are faster than "prctile" function from Statistics Toolbox from Mathworks
THR=pixels(round(length(pixels)*percentile)); 
imageCN=uint8(255*imageC/THR); % Compensated and Normalized image
%imshow(imageCN);
imageCNS=single(imageCN);
%%%


%%%%%%%%%%%% Fast template matching
frameMean = conv2(imageCNS,ones(size(kernell))./numel(kernell),'same');
templateMean = mean(kernell(:));
corrPartI = conv2(imageCNS,rot90(kernell-templateMean,2),'same')./numel(kernell);
corrPartII = frameMean.*sum(kernell(:)-templateMean);
stdFrame = sqrt(conv2(imageCNS.^2,ones(size(kernell))./numel(kernell),'same')-frameMean.^2);
stdTemplate = std(kernell(:));
corrScore = (corrPartI-corrPartII)./(stdFrame.*stdTemplate);
%%%%%%%%%%%%

%%%%%%%%%%%% Seeding and baking donuts           %%%%% SCROLL2SEED %%%%%
neuronSize = [neuronSizeIn,neuronSizeOut];
[seeds, newSeeds] = scroll2seed(corrScore, imageCN, neuronSize);
seedsigns = sign(seeds);
seedsigns(isnan(seedsigns) ) = 0;
Centers = regionprops(logical(seedsigns), 'Centroid');
%%%%%%%%%
automaticSeedsN=size(Centers,1);
for newSeedIndex=1:size(newSeeds,1)
    Centers(automaticSeedsN+newSeedIndex).Centroid=newSeeds(newSeedIndex,:);    
end
neuronR=ceil(neuronSizeOut*pixelation/2);

[maskStack,maskStackWC, maskStructure]=findDonuts(imageCN, Centers, automaticSeedsN,...
    neuronR, mode, scoreTHR,springCoeff,qualityCheckHandPicked);

imS=size(maskStack,1);
rgb(:,:,1)=imageCN;
rgb(:,:,2)=imageCN;
rgb(:,:,3)=imageCN;
rgb=double(rgb);
outerSpace=ones(size(rgb,1),size(rgb,2));

sizeoriginal = size(maskStructure);
%%%%%%%%%%%% Deselection Process [OBSOLETE]
%{      
[maskStructure, rgb, maskStack] = deselection(maskStructure, maskStack, rgb);
sizenew = size(maskStructure);
numdeleted = (sizeoriginal-sizenew);
fprintf('Deleted %d\n',numdeleted);
%}
%%%%%%%%%%
%{l
%%%%%%%%%%%%%%%%%%%%%%%% User cell selection/deselection GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% maskStructure   struct array, each index is a different cell    
% imageCN   background image matrix                              <
% outerSpace   big ass matrix                                      <
% neuronR   radius of avg Neuron?                                <
% rgbBlueCO   not defined yet                                      <
[newSeeds, maskStructure, rgb, maskStack]=playMovie( imageCN, stack, maskStructure, maskStack);
prevCenterLength = size(Centers,1);
automaticSeedsN=0;
for newSeedIndex=1:size(newSeeds,1)
    Centers(prevCenterLength+newSeedIndex).Centroid=newSeeds(newSeedIndex,:);    
end
numNewCircles = size(Centers(prevCenterLength+1:size(Centers,1)));
%%% If the user added new seeds then process them, if not do nothing
if newSeeds > 0
    % Call find donuts again using only the new coordinates added by the user
    [maskStack_N,maskStackWC_N, maskStructure_N]=findDonuts(imageCN, Centers(prevCenterLength+1:size(Centers,1)), automaticSeedsN,...
        neuronR, mode, 10 ,springCoeff,qualityCheckHandPicked);
    % We have to add the donuts generated from 'playMovie' to our current array
    % of donuts called maskStructure 
    maskStructure = [maskStructure,maskStructure_N];
    % We have to do the same for 'maskStack' so that the circles are displayed
    % properly (and also maskStackWC right after) but they are matrices so its
    % a bit longer
    A=maskStack;B=maskStack_N;
    [nRows,nCols,lenA] = size(A);
    lenB = size(B,3);
    AB = zeros(nRows,nCols,lenA+lenB);
    AB(:,:,1:lenA) = A;
    AB(:,:,lenA+1:lenA+lenB) = B;
    maskStack=AB; %
    A=maskStackWC;B=maskStackWC_N;
    [nRows,nCols,lenA] = size(A);
    lenB = size(B,3);
    AB = zeros(nRows,nCols,lenA+lenB);
    AB(:,:,1:lenA) = A;
    AB(:,:,lenA+1:lenA+lenB) = B;
    maskStackWC=AB; 
end
% Reconstructs the reference image or something I think
imS=size(maskStack,1);
rgb(:,:,1)=imageCN;
rgb(:,:,2)=imageCN;
rgb(:,:,3)=imageCN;
rgb=double(rgb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The following code adds all the blue circles to the plot
%{a
for sindex=1:size(maskStructure,2)
rgb(maskStructure(sindex).CenterCoor(1)-imS/2:maskStructure(sindex).CenterCoor(1)...
    +imS/2-1,maskStructure(sindex).CenterCoor(2)-imS/2:maskStructure(sindex).CenterCoor(2)...
    +imS/2-1,3)=rgb(maskStructure(sindex).CenterCoor(1)-imS/2:maskStructure(sindex).CenterCoor(1)...
    +imS/2-1,maskStructure(sindex).CenterCoor(2)-imS/2:maskStructure(sindex).CenterCoor(2)+imS/2-1,3) +255*maskStack(:,:,sindex);    
outerSpace(maskStructure(sindex).CenterCoor(1)-imS/2:maskStructure(sindex).CenterCoor(1)...
    +imS/2-1,maskStructure(sindex).CenterCoor(2)-imS/2:maskStructure(sindex).CenterCoor(2)...
    +imS/2-1)=outerSpace(maskStructure(sindex).CenterCoor(1)-imS/2:maskStructure(sindex).CenterCoor(1)...
    +imS/2-1,maskStructure(sindex).CenterCoor(2)-imS/2:maskStructure(sindex).CenterCoor(2)+imS/2-1) -maskStackWC(:,:,sindex);
end
%}

%figure;
imshow(uint8(rgb)) % This is the reference image with all blue selected cells
rgbBlueCO = rgb;
for y=1:length(rgbBlueCO(:,1,1))
    for x=1:length(rgbBlueCO(1,:,1))
        if rgbBlueCO(y,x,1) == rgbBlueCO(y,x,3)
            %rgbBlueCO(y,x,1)=1;
            %rgbBlueCO(y,x,2)=1;
            %rgbBlueCO(y,x,3)=1;
        else
            %rgbBlueCO(y,x,1)=1;
            %rgbBlueCO(y,x,2)=1;
            rgbBlueCO(y,x,3)=rgbBlueCO(y,x,3)/1.8;
        end
    end
end
% pause;
imshow(uint8(rgb))
% imshow(outerSpace);
% pause;

saveas(gca,[fileName,' clear reference image.png']);

end