function [  ] = intensityZ( stackAdjustedCut, intavgNPCorr, firingTimes,...
                            binaryFiring, firedNeurons, file)
%INTENSITYZ Summary of this function goes here
%   Detailed explanation goes here

outputName = [file,' final'];

fps = 3.91;
frames = size(intavgNPCorr,2); 
totalCells = size(intavgNPCorr,1);
activeCells = length(firedNeurons);
seconds=round(frames/fps);
lag=10;

%%%%%%%%% Creates a Raster Plot of the data %%%
myRasterPlot(binaryFiring, firedNeurons, file);

%%%%%%%%% Plots a histogram of cells firing over time %%%
hold off
hist(1./firingTimes(firingTimes<1),frames)
title(['Digitized events over time, ',num2str(totalCells), ' total cells.'])
xlabel('Frame number')
ylabel('Number of cells firing')
pause(1.5);
saveas(gcf,[outputName,' Digitized cell events.png'])

%%%%%%%%% Plots percent of cells fired over time %%%
binSize = 4; % bin size in frames, about 1 second.
numActive = zeros(ceil(frames/binSize)); % Crunches frames down to arbitrary step size
percentActive = zeros(ceil(frames/binSize));
fireTimesSingleEvent = sort(round(1./firingTimes(firingTimes<1)));
for i=1:frames
    numActive(ceil(i/binSize)) = sum(fireTimesSingleEvent==i) + numActive(ceil(i/binSize));
 %   percentActive(ceil(i/binSize)) = percentActive(ceil(i/binSize)) + numActive(i);
end
percentActive = 100*numActive/activeCells;
hold off
plot(percentActive);
title(['Percent average cells over time, bin size = ',num2str(binSize),' frames'])
xlabel('Bins (Seconds if bin size is 4)')
ylabel('Percent of active cells firing')
pause(1.5);
saveas(gcf,[outputName,' percent active.png'])

%%%%%%%%% Plots intensity over time %%%
hold off;
intensityCellFrames = zeros(1,frames);
intensityEntireImage = zeros(1,frames);
% Calculates the average intensity of the cells and the whole image for each frame
for i=1:frames
   intensityCellFrames(i) = mean(intavgNPCorr(:,i));
%    intensityEntireImage(i) = mean(mean(stackAdjustedCut(:,:,i)));
end
%%% Applies a exponential fit to the intensity data, doesn't use curve fitting toolbox
fit = @(b,x) b(1).*exp(-.01*(b(2)).*(x.^2));
x  = linspace( 1, frames, frames); y = intensityCellFrames;
OLS = @(b) sum((fit(b,x) - y).^2);          % Ordinary Least Squares cost function
opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
B = fminsearch(OLS, rand(2,1), opts);       % Use ?fminsearch? to minimise the ?OLS? function

avgIntensity = mean(intensityCellFrames);
for c=1:frames
    %intensityData(r,c) = (intensityData(r,c)-fit(B,c))/fit(B,c); %%% Uncomment this for Exp fit
%    intensityCellFrames(c) = (intensityCellFrames(c)-fit(B,c))/fit(B,c);
end
intensityCellFrames = intensityCellFrames-avgIntensity;
T = 1; t = 1/0.02; a = T/t;
% where a = T/?, T = the time between samples, and ? (tau) is the filter time constant.
xfilt = filter([1-a a-1],[1 a-1], intensityCellFrames)+500;

% plot(intensityCellFrames,'r')
plot(xfilt,'r','DisplayName','Filtered intensity (HP)')
hold on
% Plot a running average of the raw data & the raw data
laggedIntensityCellFrames = conv(intensityCellFrames, ones(1,lag), 'valid')/lag;
plot(laggedIntensityCellFrames,'k','linewidth',2,'DisplayName','Lagged intensity over time')
plot(intensityCellFrames,'b','linewidth',.1,'DisplayName','Intensity over time')
legend('show')
% This line will plot a running average for xfilt
%laggedxfilt = conv(xfilt, ones(1,lag), 'valid')/lag; plot(laggedxfilt,'k','linewidth',2)

stdDev = std(intensityCellFrames);
stdDevFilt = std(xfilt);
%plot(intensityEntireImage)
title(['Average fluorescent intensity across all cells over time, stdDev=', num2str(stdDevFilt)])
xlabel('Frame number')
ylabel('Change in F over polyfit baseline F')
saveas(gcf,[outputName,' Deviation from average intensity per frame.png'])



%%%%%%%%% Histogram Plot of total deltaF/F fluctuations


end
