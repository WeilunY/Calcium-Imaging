function [ ] = motionCutterGUI( motionCompensation, stackAdjusted,Y1,outputName,motionData)
%MOTIONREMOVE Summary of this function goes here

%motionData = true; %%% Set to true if you DO have a .mat file with motion data

slmin=1;slmax=40;
sliderVal=20;
threshold_factor=sliderVal;
uicontrol('Style','slider','Callback',@sliderCallback,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderVal,...
                    'Position',[5 5 200 20]);
init_ui();

function sliderCallback(src, evt)
    sliderVal=round(get(src, 'Value'));
    threshold_factor = sliderVal;
    init_ui();
end

function init_ui()
    uicontrol('Style','slider','Callback',@sliderCallback,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderVal,...
                    'Position',[5 5 200 20]);
    % lag indicates the number of previous data points used with the
    % current data point when calculating the moving average
    lag = 5;
    % minCut is the minimum number of frames that can be cut together,
    % this is so we cut in clusters not single frames
    minCut = 5;% * NOT YET IMPLEMENTED, PROBABLY NOT NECESSARY *

    originalCorr = motionCompensation(1,:);% Correlation of original images to ref image
    correctedCorr = motionCompensation(2,:);% Correlation of corrected images to ref image
    hold off
    %plot( originalCorr,'r','linewidth',.5 )
    plot( correctedCorr )
    hold on
    
    % Filter stuff
%     input = correctedCorr;
%     d = fdesign.lowpass('Fp,Fst,Ap,Ast',3,5,0.5,40,4);
%     Hd = design(d,'equiripple');
%     output = filter(Hd,input);
%     fvtool(Hd)
%     plot(output,'r')
    % End filter stuff
    
    % Calculate avg correlation and stdDev of correlation, use them to
    %   decide on the threshold
    avgCorrCoeff = mean(correctedCorr);
    stdDevCorrCoeff = std(correctedCorr);
    % Any frames with an avg correlation below this will be cut
    threshold = (avgCorrCoeff-stdDevCorrCoeff)+(threshold_factor-20)/100
    
    % Code for plotting a horozontal line (threshold value)
    L_min = 0;
    L_max = size(correctedCorr, 2);
    y = threshold; %constant y value
    d = 1; %direction
    r = 15; %y range
    a = linspace(L_min, L_max, 300);
    b = linspace(y, y, 300);
    plot(a,b,'--k','linewidth',1.5)
    % Done plotting horozontal line
    
    % CHANGE: Add current threshold value next to the line
    txt = ['Current Threshold: ' num2str(y)];
    text(400, y + 0.02, txt);
    

    % Moving avg for original image
 %   avgOriginalCorr = conv(originalCorr, ones(1,lag), 'valid')/lag; (not used)
    % Moving avg for corrected image
    avgCorrectedCorr = conv(correctedCorr, ones(1,lag), 'valid')/lag;

    %plot(avgOriginalCorr,'--')
    plot(avgCorrectedCorr,'k','linewidth',2)
    title('Moving average of motion correlation and threshold')
    xlabel('Frames')
    ylabel('Correlation Coefficient')
     hold off
    %saveas(gcf,[outputName,' motion cut threshold.png'])
    %pause(1.5);
end


pause();
saveas(gcf,[outputName,' motion cut threshold.png'])
close();

% We go through the correlation array and record frames with low correlation
framesCut=find(avgCorrectedCorr<threshold);
framesNotCut=find(avgCorrectedCorr>threshold);

stackAdjustedCut=stackAdjusted;%%% The new cut stack
stackAdjustedCut(:,:,framesCut)=[];

cutfromStack=stackAdjusted;%%% all frames that were cut
cutfromStack(:,:,framesNotCut)=[];

frames = size(avgCorrectedCorr, 2);
frameCorrStatus = ones(1,frames);
frameCorrStatus(framesCut)=0;

%%% This chunk will make a plot of the mouse's movement over time, removing
%%% the sections where it moves too much
if motionData
Y11=Y1(:,1);  
compareCorr = zeros(length(Y11),1);
lengthRatio = floor(length(Y11)/length(frameCorrStatus));
%%%%% Bullshit line makes this work somehow  lengthRatio = lengthRatio - 250;
for i=1:length(frameCorrStatus)-2%%% Must resize frame Corr Status to Y11 dimensions
    for q=1:lengthRatio
        compareCorr(lengthRatio*i+q)=frameCorrStatus(i);
        if frameCorrStatus(i)==1 && (lengthRatio*i+q)<size(Y11,1)%%% If frame isn't cut
            %Y11(lengthRatio*i+q)=-1;
            Y11(lengthRatio*i+q)=Y11(lengthRatio*i+q)-125; % ERROR for some reason??
        end
    end
end

%%% Now I plot the data and save it
%{
plot(Y1(:,1))      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
plot(Y11)
plot(Y11(compareCorr > 0)-200)
legend('Uncut movement','Overlayed part has been CUT','New, cut movement')
% This chunk will create a second set of Axis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax1_pos = get(gca,'Position'); % position of first axes
ax2 = axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
% sets the y tick labels
set(ax2,'YTick',0:1,'YTickLabel',['' ''])
%seconds=round(length(Intensity_DeltaF(1,:))/3.91);
% sets the x tick labels
set(ax2,'XTick',0:2,'XTickLabel',[0 length(motionCompensation)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
title('Plot of mouse movement before and after cutting')
%}
%saveas(gcf,[outputName,' cut plot of movement.png'])
Y11cut=Y11(compareCorr > 0);
close all
hold off
save(outputName, 'frameCorrStatus', 'stackAdjustedCut','cutfromStack','Y11cut')
end
%pause();

save(outputName, 'frameCorrStatus', 'stackAdjustedCut','cutfromStack')
load([outputName,'.mat']);
end



% motionCutter( motionCompensation, stackAdjusted,'nem' )

