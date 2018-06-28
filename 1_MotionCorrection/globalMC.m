function [stackOriginal, stack, stackAdjusted, motionCompensation, timeAxis, normalizedMotion]=...
    globalMC(FileName1,PathName1,FileName2,PathName2,motionData)
%%%%%%%%%%%%%%%%%% Parameters
%[FileName1,PathName1] = uigetfile({'*.tif'},'Select your time series data')%%%Data source
%[FileName2,PathName2] = uigetfile({'*.mat'},'Select the motion data exported from WCP')%%%Data source
%%%%%%%%%%%%%%%%%%%%%% Calm period parameters
%motionData = true; %%% Set to true if you DO have a .mat file with motion data

if motionData
    load([PathName2, FileName2]);
    movementSF=round(length(T1)/T1(length(T1))); % Sampling rate of the WCP record
else
    movementSF=112.6; % Assumed
    fps=3.91;
end

% This is a parameter of joining range. If there is a problem with finding a calm period, increase this parameter
imdilateParameter=12; 
imdilateParameter=12;

%%%%%%%%%%%%%%%%%%%%%%Reference finder parameters
Nimages=100; %%%% How many images from the calm period should be each2each correlate
Corrperc=0.95; %%%% What part of the distribution of slide2slide corelation 
                 %%%  should be taken to pick the best aligned frames from the calm period. 
%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%% Registration parameters
dispDynamicCoeff=0.5;
dispHardCoeff=1-dispDynamicCoeff;
maxDisplacement=8;% maximal displacement (from the reference) movements. In pixels;
penaltyPow=0.01; % Power of the interrow displacement penalty. The higher value the harder to slide subsequent rows on each other.
%It is further multiplied by a "normalizing" factor derived from the
%average inter line correlation within a single image
%%%%%%%%%%%%%%%%%%%%%%

[stack]=reader([PathName1 FileName1]); %%%%%% Image loading
imagesN=size(stack,3);

%%%%%%%%%%%%%%%%%%Find reference for registration
%%% If motion data exists, find calm period based on it.
if motionData 
    [calmTimeSt,calmTimeEnd,expStarted,expEnded, Y1, M_calm, M_mov]=...
        calmPeriod(movementSF, [PathName2 FileName2],imdilateParameter);
    M_calm=M_calm-expStarted;%%%%%%%% Bounds for dynamic displacement Limit
    M_mov=M_mov-expStarted;
else
    expStarted = 0;  
    expEnded = 299.9; % Assuming 5 minute file
%    calmTimeSt = 315; %%%% Calm Period Start in Seconds <- changes with file
%    calmTimeEnd = 460;  %%%% Calm Period End in Seconds   <- changes with file
%   M_calm=[2 2 2 2 2 2 2 2 2 2];
%   M_mov=[2 2 2 2 2 2 2 2 2 2];
    [calmTimeSt,calmTimeEnd, M_calm, M_mov] = calmPeriodNoMD(stack);
    M_calm=M_calm-expStarted;%%%%%%%% Bounds for dynamic displacement Limit
    M_mov=M_mov-expStarted;
end
fps=3.91;% imagesN/(expEnded-expStarted);  %%% Full frame sampling rate;

[reference]=referenceFinder(stack, fps, expStarted, calmTimeSt, calmTimeEnd, Nimages, Corrperc);

%%%%%%%%%%%%%%%Corection of the reference for ultrabright neurons
brightestPerc=0.95; %%%Brightness level to be clipped

referenceCropped=reference;
pixels=sort(referenceCropped(:));
THRgray=pixels(round(length(pixels)*brightestPerc));
body=find(reference>THRgray);
referenceCropped(body)=THRgray;

%%%%%%%% Bounds for dynamic displacementLimit
M_max=zeros(1,10);
M_min=zeros(1,10);

size(stack);
for index=1:10
    if (M_calm(index)-1)*fps <= 0
        M_calm(index) = M_calm(index-1) + 1;
    end
    rawImage=stack(:,:,round((M_calm(index)-1)*fps));
    [~, ~, M_max(index)]=globalShifter2(rawImage, referenceCropped, maxDisplacement); 
    rawImage=stack(:,:,round(M_mov(index)*fps));
    [~, ~, M_min(index)]=globalShifter2(rawImage, referenceCropped, maxDisplacement);  
end

M_max=mean(M_max);
M_min=mean(M_min);
M_range=M_max-M_min;
%%%%%%%%
%%%%%%%%%%%%%%%Image registration
stackAdjusted=zeros(size(stack), 'uint16');

h = waitbar(0,'Correcting for Motion Frame by Frame');
motionCompensation1=zeros(1,imagesN);
motionCompensation2=zeros(1,imagesN);
% This is bypass for dynamic maxDisplacement. Choose as you wish.
maxDisplacementHARD=maxDisplacement;
tic

% This for loop traverses the tiff stack and corrects each frame
for index=1:imagesN
    waitbar(index/ imagesN)
    rawImage=stack(:,:,index);
    
    %%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%% removes extremely bright parts of 'rawImage'
    %%% This has a bug, round(length(pixels)*brightestPerc) is not correct,
    %%% it should be the value of the brightest pixel
%   rawImageCropped=rawImage;
%   pixels=sort(rawImageCropped(:));
%   THRgray=pixels(round(length(pixels)*brightestPerc));
%   points=find(rawImage>THRgray);
%   rawImageCropped(points)=THRgray;
%   rawImage = rawImageCropped;
    %%% This code drastically reduces image quality, can't figure out why %%%
    %%%%%%%%%%%%%%%%%%%%%%% TEST CODE %%%
    
    %A "normalizing" factor derived from the average inter line correlation within a single image
   rByRcorr1=mean(diag(corrcoef(double(rawImage)),1)); 
   rByRcorr2=mean(diag(corrcoef(double(rawImage')),1));
    
    [~, adjustedImage, ~]=globalShifter2(rawImage, referenceCropped,maxDisplacement);
    
%maxDisplacementTMP=abs(ceil(maxDisplacementHARD*(dispHardCoeff+dispDynamicCoeff*((M_max-M)/M_range))));
    [~, adjustedImage]=rowShifter2(adjustedImage, referenceCropped, maxDisplacementHARD, abs(rByRcorr1), penaltyPow);
    [~, adjustedImage]=rowShifter2(adjustedImage', referenceCropped',maxDisplacementHARD, abs(rByRcorr2), penaltyPow);

    adjustedImage=adjustedImage';
    stackAdjusted(:,:,index)=adjustedImage;

    % Correlation of original images to reference image.
motionCompensation1(1,index)=corr2(rawImage(2*maxDisplacement+1:end-2*maxDisplacement,2*maxDisplacement+1:end-2*maxDisplacement),...
    referenceCropped(2*maxDisplacement+1:end-2*maxDisplacement,2*maxDisplacement+1:end-2*maxDisplacement)); 
    % Correlation of corrected images to reference image.
motionCompensation2(1,index)=corr2(adjustedImage(2*maxDisplacement+1:end-2*maxDisplacement,2*maxDisplacement+1:end-2*maxDisplacement),...
    referenceCropped(2*maxDisplacement+1:end-2*maxDisplacement,2*maxDisplacement+1:end-2*maxDisplacement));


end

toc
close(h);
motionCompensation=[motionCompensation1; motionCompensation2];
stackOriginal=stack;
stack(end-maxDisplacement+1:end,:, :)=[];
stack(1:maxDisplacement,:, :)=[];
stack(:,end-maxDisplacement+1:end, :)=[];
stack(:,1:maxDisplacement, :)=[];
stackAdjusted(end-maxDisplacement+1:end,:, :)=[];
stackAdjusted(1:maxDisplacement,:, :)=[];
stackAdjusted(:,end-maxDisplacement+1:end, :)=[];
stackAdjusted(:,1:maxDisplacement, :)=[];
a = waitbar(1,'Saving the results. Large stacks.');

if motionData
    normalizedMotion=Y1(round(linspace(expStarted*movementSF, expEnded*movementSF,imagesN))+1)/mean(Y1(:,1));
else
    normalizedMotion=0;
end

timeAxis=linspace(expStarted,expEnded, imagesN);
FileName1(end-2:end)='mat';
fpath = [PathName1, 'Data/'];
save(fullfile(fpath, FileName1),'stackOriginal', 'stack', 'stackAdjusted', 'motionCompensation', 'timeAxis', 'normalizedMotion');
close (a);

plot(timeAxis,normalizedMotion,'DisplayName','normalizedMotion');
hold on
plot(timeAxis,motionCompensation(1,:),'r','DisplayName','Correlation of original images to reference image');
plot(timeAxis,motionCompensation(2,:),'g','DisplayName','Correlation of corrected images to reference image');
legend('show')

% Correlation of original images to reference image.
% Correlation of corrected images to reference image.
grid minor
FileName1(end-2:end)='fig';
fpath = [PathName1, 'Data/'];
saveas(gcf, fullfile(fpath, FileName1), 'fig');
%pause;
close;
end