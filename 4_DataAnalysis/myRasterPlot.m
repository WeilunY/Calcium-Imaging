function [ ] = myRasterPlot( binaryFiring, firedNeurons, fileName )
% This creates a raster plot of inputted data

numPlots = length(firedNeurons);
frames = size(binaryFiring,2);
seconds=round(frames/3.91);

% plotting the actual data
hold off; % clear remaining figures
clf
yPos = 1;
% Cell i
for i=1:numPlots
    % Frame ii
    for ii=1:frames
        if binaryFiring(i,ii)==1
            plot([ii,ii+1], [yPos yPos],'--k','linewidth',2.7);
        end
    end
%    plot( binaryFiring(i,:)+(inc-1),'linewidth',1.1)
    yPos = yPos+1;
    hold on;
end
% done plotting data

fontsize = 8;
if numPlots > 150
    fontsize = 5;
end

axis([0,frames,-0.5,numPlots-0.2])
title(['\fontsize{12}Raster Plot of ',strrep(fileName,'_',' ')])
xlabel('\fontsize{12}Seconds')
ylabel('\fontsize{12}Cell Number')

%%% Setting the x tick marks
xticks([0 frames*0.1 frames*0.2 frames*0.3 frames*0.4 frames*0.5 ...
    frames*0.6 frames*0.7 frames*0.8 frames*0.9 frames])
xticklabels({0 seconds*0.1 seconds*0.2 seconds*0.3 seconds*0.4 seconds*0.5 ...
    seconds*0.6 seconds*0.7 seconds*0.8 seconds*0.9 seconds})

%%% Setting y tick marks, will remove some if too many plots
N = round(numPlots/100);
firedNeuronsCut = firedNeurons(1:N:length(firedNeurons));
set(gca,'YTick',0:N:numPlots-1,'YTickLabel',firedNeuronsCut(1:length(firedNeuronsCut)),'FontSize',fontsize)
%set(gca,'XTick',tickVals,'XTickLabel',tickText)

saveas(gcf,[fileName,' Raster Plot'])
end

