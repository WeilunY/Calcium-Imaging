function [maskStack, maskStackWC, maskStructure]=findDonuts(imageCN,...
    Centers, automaticSeedsN, neuronR, mode, scoreTHR, springCoeff,qualityCheckHandPicked)
%imageCN ... compensated and normalized image
%Centers ... centers of identified neurons
%neuronR ... readius of neuron in pixels
%mode ... type of neuronal shape and boundaries - =1 for forced circular
%neurons (faster, robust), =2 for irregular neurons with some constraints
%(slower, possible bugs, but theoretically most accurate), 3= irregular
%neurons and irregular boundaries
% automaticSeedsN is The number of automatically selected cells are in
%     'Centers', the higher indices were selected by hand1

minDist=20; % The minimum seperation distance between the centers of the circles (in pixels)
sigma=max(round((0.125)*neuronR),1);

H = fspecial('gaussian',3*sigma,sigma);
imageCN_blurred=imfilter(imageCN,H);
neuronField=round(neuronR*1.333);
neuronField2=2*neuronField;

%%% Polar  - Cartesian transformations
[theta,rho]=meshgrid(linspace(0,(71/36)*pi,72), neuronField:-1:1); % 5° step in angle
[X,Y] = pol2cart(theta,rho);
xx=round(X);
yy=round(Y);
linearInd = sub2ind([2*neuronField+1 2*neuronField+1], xx(:)+neuronField+1, yy(:)+neuronField+1);
%%%
[rN, cN]=size(imageCN);

%%%% Treatment of cutted neurons (TODO Spare the poor boundary guys)
cuttedNeurons=[];
if qualityCheckHandPicked==0
    for indexCutted=1:automaticSeedsN
        xcor = Centers(indexCutted,1).Centroid(1);
        ycor = Centers(indexCutted,1).Centroid(2);
        if round(xcor+neuronField)>cN || round(xcor-neuronField)<1 ||...
                round(ycor+neuronField)>rN || round(ycor-neuronField)<1
                cuttedNeurons=[cuttedNeurons indexCutted];
        else
            for nestedIndex=indexCutted+1:size(Centers,1)
                xcorNew = Centers(nestedIndex,1).Centroid(1);
                ycorNew = Centers(nestedIndex,1).Centroid(2);
                if (xcorNew-xcor)^2+(ycorNew-ycor)^2 < minDist^2
                    cuttedNeurons=[cuttedNeurons indexCutted];
                end
            end
        end
    end
else
    for indexCutted=1:size(Centers,1)
        xcor = Centers(indexCutted,1).Centroid(1);
        ycor = Centers(indexCutted,1).Centroid(2);
        if round(xcor+neuronField)>cN || round(xcor-neuronField)<1 ||...
                round(ycor+neuronField)>rN || round(ycor-neuronField)<1
            cuttedNeurons=[cuttedNeurons indexCutted];
        else
            for nestedIndex=indexCutted+1:size(Centers,1)
                xcorNew = Centers(nestedIndex,1).Centroid(1);
                ycorNew = Centers(nestedIndex,1).Centroid(2);
                if (xcorNew-xcor)^2+(ycorNew-ycor)^2 < minDist^2
                    cuttedNeurons=[cuttedNeurons indexCutted];
                end
            end
        end
    end
end
autoNdeleted=numel(find(cuttedNeurons<=automaticSeedsN));
Centers(cuttedNeurons)=[];
innerNeuronsN=size(Centers,1);

%%%%
maskStack=zeros(neuronField2,neuronField2,innerNeuronsN);
maskStackWC=zeros(neuronField2,neuronField2,innerNeuronsN);
%%%%
phi=linspace(0,pi*71/36,72);
CoorX=cos(phi);
CoorY=sin(phi);

scoreRecord=zeros(innerNeuronsN,1);

refinedIndex=1;
maskStructure(innerNeuronsN).CenterCoor=[];
maskStructure(innerNeuronsN).surroundSize=[];
maskStructure(innerNeuronsN).WCPixList=[];
maskStructure(innerNeuronsN).nucleusPixList=[];
maskStructure(innerNeuronsN).ringPixList=[];
autoRemoved = 0;
minradius = 9999;
maxradius = 0;
for index=1:innerNeuronsN
    currNeuronIm=imageCN_blurred(max(round(Centers(index,1).Centroid(2))...
        -neuronField,1):min(round(Centers(index,1).Centroid(2))+...
        neuronField,rN),max(round(Centers(index,1).Centroid(1))-...
        neuronField,1):min(round(Centers(index,1).Centroid(1))+neuronField,cN));
    %%%%%k èemu jsou tady ty min max???TODO zbyteèné?
    %%%%% I think I was testing something here with the next two lines? Maybe... Keep eye on it
%     size(currNeuronIm)
%     linearInd
    metaRepre=currNeuronIm(linearInd);
    currNeuronImPolar=reshape(metaRepre, [neuronField 72]);
    %[Centers(index,1).Centroid(2) Centers(index,1).Centroid(1)]
    rProjection=mean(currNeuronImPolar,2);
    ring=mean(currNeuronImPolar(round(neuronField/4):round(neuronField/2),:),1);
    %%% Score and clear false positive neurons
    normFactor=mean2(currNeuronImPolar(round(neuronField/4):neuronField,:));
    %neuronR
    ringMean=mean(ring);
    ringStd=std(ring);
    rProjectionMean=mean(rProjection);
    
    ringScore=(ringMean-rProjectionMean)/normFactor;
    nucleusScore=(ringMean-mean(rProjection(neuronR:end)))/normFactor;
    outerScore=(ringMean-mean(rProjection(1:round((1/3)*neuronR))))/normFactor;
    smoothRingPenatly=ringStd/ringMean;
    
    %!!!!!!!!!!!!Score function
    %%%This is the neuron-vs.-splodge evaluation function. Please modify, or supply with coefficients if necessary.
    scoreFcn=ringScore+nucleusScore+outerScore-smoothRingPenatly; 
    scoreRecord(index)=scoreFcn;
    

    if (scoreFcn>scoreTHR)||(index>(automaticSeedsN-autoNdeleted))
    %%% Visualisation of a suspected area
%     subplot(2,2,1);
%     imshow(currNeuronImPolar,'InitialMagnification',400);
%     subplot(2,2,2);
%     plot(rProjection/max(mean(currNeuronImPolar,2)));
%     ylim([0 1])
% grid minor;
%     subplot(2,2,3);
%     imshow(currNeuronIm);
%     subplot(2,2,4);
%     plot(ring);
%     ylim([0 255])
%     grid minor;
        
        %%% According to mode callculate and write down the donuts
        switch mode
            case {1}
                %TODO circular donuts with a constant thicknes
            case {2}
                %Irregular neurons with a constant thicknes
                [ridgeRows]=ridgeFinder(currNeuronImPolar,springCoeff);
                %%% Added by Alex %%% Removes any cell not within a certain radius threshold
                avgRadius = 150/mean(ridgeRows); % 'radius' in arbitrary weight units
                minradius = min([minradius,avgRadius]);
                maxradius = max([maxradius,avgRadius]);
                %if avgRadius < 13 || avgRadius > 20
                if avgRadius < 12
                    continue    %%% Ignores very small circles
                end
                innerOuterEdges=[0.8;1.2]*(neuronField-ridgeRows+1);
                innerOuterCoorCol=round(innerOuterEdges(1,:).*CoorX);
                innerOuterCoorCol(2,:)=round(innerOuterEdges(2,:).*CoorX);                
                innerOuterCoorRow=round(innerOuterEdges(1,:).*CoorY);
                innerOuterCoorRow(2,:)=round(innerOuterEdges(2,:).*CoorY);
                nucleus  = poly2mask(innerOuterCoorCol(1,:)+neuronField, ...
                    innerOuterCoorRow(1,:)+neuronField, neuronField2, neuronField2);
                wholeCell = poly2mask(innerOuterCoorCol(2,:)+neuronField, ...
                    innerOuterCoorRow(2,:)+neuronField, neuronField2, neuronField2);
%                  subplot(2,2,1);
%                 imagesc(nucleus);
%                 subplot(2,2,2);
%                 imagesc(wholeCell);

                maskStack(:,:,refinedIndex)=wholeCell-nucleus;
                maskStackWC(:,:,refinedIndex)=wholeCell;
%                 wholeCellInd=find(wholeCell);
%                 nucleusInd=find(nucleus);
                    %close all;
            case {3}
                %TODO irregular neurons with a threshold-based thicknes
        end
        maskStructure(refinedIndex).CenterCoor=[Centers(index,1).Centroid(2) Centers(index,1).Centroid(1)];
        maskStructure(refinedIndex).surroundSize=size(currNeuronIm);
        maskStructure(refinedIndex).WCPixList=find(wholeCell);
        maskStructure(refinedIndex).nucleusPixList=find(nucleus);
        maskStructure(refinedIndex).ringPixList=setdiff(maskStructure(refinedIndex).WCPixList,...
            maskStructure(refinedIndex).nucleusPixList);
        refinedIndex=refinedIndex+1;
    end    

end
minradius;
maxradius;

maskStack(:,:,refinedIndex:end)=[];
maskStackWC(:,:,refinedIndex:end)=[];
maskStructure(refinedIndex:end)=[];
% hist(scoreRecord,40);
% pause

end