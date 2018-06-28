function [ ] = intensityFullPlot( Intensity_DeltaF,firedNeurons,firingThresholds )
%INTENSITY_DELTAF_PLOTTER plots deltaF intensityData

% intensityFullPlot( intensityData,firedNeurons,firingThresholds )

numPlots=length(firedNeurons);
remainingCells=[];

% clear remaining figures
clf
for i=1:(numPlots)
    plot( Intensity_DeltaF(firedNeurons(i),:)+(i-1),'linewidth',1.1)
    hold on;
end
% done plotting data

%code for plotting a horozontal line
L_min = 0;
L_max = length(Intensity_DeltaF(1,:));
%L_min = L_max/12;
y = 0; % constant y value
d = 1; % direction
r = 15;% y range
a = linspace(L_min, L_max, 2);
b = linspace(y, y, 2);
plot(a,b,'--k','linewidth',1.5)
for i=1:((numPlots)*2-1)
    y = i/2;
    if mod(i,2)==1%%% If this horozontal line is a threshold
        y = round(y-0.5)+firingThresholds(round(i/2));
        a = linspace(L_min, L_max, 2);
        b = linspace(y, y, 2);
        plot(a,b,'--k','linewidth',1.2)
    else%%% If line is an avg
        a = linspace(L_min, L_max, 2);
        b = linspace(y, y, 2);
        plot(a,b,'--k','linewidth',0.8)
    end
end
% done plotting horozontal lines

% Should eliminate white space from the plot (only works w/1 plot)
%set(gca,'LooseInset',get(gca,'TightInset'))
% Chooses range of the plot
axis([0,L_max,-0.5,numPlots-0.5])
% Sets the y tick labels
set(gca,'YTick',0:numPlots-1,'YTickLabel',firedNeurons)
title('Cell Intensities with Firing Thresholds')
ylabel('Cell Number')
xlabel('Frames')

end