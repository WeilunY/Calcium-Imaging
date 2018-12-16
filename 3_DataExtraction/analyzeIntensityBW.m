function [  ] = analyzeIntensityBW(intavgNPCorr, outputName)
ignoreInitialSpiking = false; % Line 140 to change analyzable part
stdCutoff = 0.053;
stdCutoff = 0.065;

%%% DeltaF/F %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%
%{
intavgNPCorrT = intavgNPCorr;
for cell=1:initialTotalCells
    for frame=2:length(intavgNPCorr(1,:))
        currVal = intavgNPCorr(cell,frame);
        prevVal = intavgNPCorr(cell,frame-1);
        intavgNPCorrT(cell,frame) = (currVal-prevVal)/currVal;
    end
    intavgNPCorrT(cell,1) = intavgNPCorrT(cell,2);
end
intavgNPCorr = intavgNPCorrT;
%plot (intavgNPCorr(1,:))
%pause()
%}
%%% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %% %%

%%% Removes "cells" with an average intensity below a threshold
%%% indicating that it might be just dark neuropil
darknessThreshold = 900; % Above 1000 seems to indicate a real cell
darknessThreshold = 100;
r=1;
while r < length(intavgNPCorr(:,1))
    darkLevel = mean(intavgNPCorr(r,:));
    if darkLevel < darknessThreshold
        intavgNPCorr(r,:)=[];
        r=r-1;
        fprintf('Removing a cell, too dark, light level: %d\n',darkLevel);
    end
    r=r+1;
end

fps = 3.91;
frames = length(intavgNPCorr(1,:));
totalCells = length(intavgNPCorr(:,1));
totalCells
realCells = length(intavgNPCorr(:,1));
lag = 31; %%% If lag is 0 then a moving average is NOT used
intensityData = intavgNPCorr;
movingavgIntensity=[];
% Multiplied by negative standard deviation is the threshold for firing
threshold_factor = 8.0; % threshold = std_neg * threshold_factor;

x = linspace(1,frames,frames);
for r=1:totalCells
    y = intensityData(r,:);
    
    %%% Exponential fitting function [NO CURVE FITTING LIBRARY]
%     fit = @(b,x) b(1).*exp(-b(2).*x);             % Objective function
%     fit = @(b,x) b(1).*exp(-.01*(b(2)).*(x.^2));
%     x  = linspace( 1, frames, frames);
%     OLS = @(b) sum((fit(b,x) - y).^2);          % Ordinary Least Squares cost function
%     opts = optimset('MaxFunEvals',50000, 'MaxIter',10000);
%     B = fminsearch(OLS, rand(2,1), opts);       % Use ?fminsearch? to minimise the ?OLS? function
%     for c=1:length(intensityData(1,:))
%         intensityData(r,c) = (intensityData(r,c)-fit(B,c))/fit(B,c);
%     end
	%%% Exponential fitting function [Uses library, better]
    x = (1:1:size(y'));
    f = fit(x',y','exp2','Start',[1000,-0.001,3500,0.0001]);
    for c=1:length(intensityData(1,:))
         intensityData(r,c) = (intensityData(r,c)-f(c))/f(c);
    end
    
%     if r==totalCells
        %hold on
        %plot(fit(B,x))
        %plot(intensityData(r-1,:))
        %pause()
%     end
end

assignin('base', 'Intensity_DeltaF', intensityData);
firedNeurons=[]; % An array that holds each firing cell's number
firingThresholds=[]; % The minimum threshold the signal must reach
firingTimes=[]; % Holds the frame of every cell that fires
%%% firingTimes is weird. It holds both the cell numbers and their firing times.
%%% #>1 is cell number, followed by (1/frameFired). All numbers greater
%%% than one is a cell number and every number less than one is the
%%% inverted frame on which it fired
binaryFiring = zeros(totalCells,frames);
cellPlot = 1;
numCellPlots = 1;

f = figure('KeyPressFcn',@keypress);
Calculate_Events();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
slmin=4.5;slmax=10;
uicontrol('Style','slider','Callback',@sliderCallback,'Min',slmin,'Max',slmax,...
                    'SliderStep',[.5 .5]./(slmax-slmin),'Value',threshold_factor,...
                    'Position',[5 5 200 20]);
sliderVal=threshold_factor;

function keypress(src, evt)
    val = double(get(f,'CurrentCharacter'));
    switch val%%%% REPLACE THIS CRAP WITH ANOTHER SLIDER THAT CHOOSES WHICH PLOT
        case 28 % <-
            if cellPlot <= 1
                fprintf('Already on first plot');
                cellPlot=1;
            else
                cellPlot=cellPlot-1;
                Intensity_DeltaF_Plotter();
            end
            init_ui();
        case 29 % ->
            if cellPlot >= numCellPlots
                fprintf('Already on last plot');
            else
                cellPlot=cellPlot+1;
                Intensity_DeltaF_Plotter();
            end
            init_ui();
        otherwise
            fprintf('Button not recognized, value is: %d\n',val);
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = Calculate_Events()
    cellPlot=1;
    std_neg = 0;
    amplitudes = [];
    std_negVals = ones(totalCells,1);
    
    number_of_events = 0;
    firingThresholds = [];
    firingTimes=[];
    avgEventSizes = [];
    
%%% For every segmented cell (cell i) run intensity analysis
%%% cell i
    for i=1:length(intensityData(:,1))%%% Data analysis
        eventsize=0;
        cellFireTimes = (i); % Keeps track of the times at which the cell fires
        A = intensityData(i,:);
        std_neg = std(A(A<0)); % Negative Std Dev of the signal
        std_negVals(i) = std_neg;
        threshold = std_neg*threshold_factor;
        event_vals=[]; % An array holding every single event size
        
        startFrame = 5;
        if ignoreInitialSpiking
            startFrame = round(length(A)/16);
        end
        if std_neg>stdCutoff
            realCells=realCells-1;
            realCells;
        end
    %%% Registers an event only if a value above the threshold is preceded by a value below it
    %%% frame ii
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ii=startFrame:(length(A)-2) 
    %%% The following if statement will determine if a neuron has fired
            if (A(ii)<threshold) && (A(ii+1)>threshold) && std_neg<stdCutoff %&& (A(ii+2)>0)
                number_of_events=number_of_events+1;
                cellFireTimes = [cellFireTimes,1/(ii+1)];
                amplitudes(number_of_events) = A(ii+1);
                iii=ii+1;
    %%% Calculating size of the event, add until it is below the threshold
    %%% subframe iii
                while (A(iii)>threshold/3 && iii<length(A))
                    event_vals=[event_vals, A(iii)];
                    binaryFiring(i,iii)=1;
                    iii=iii+1;
                end
                %%% 800ms x 3.91fps =~= 3 frames, events cannot be closer
         %      ii=ii+round(.8*fps);
             
                ii=ii+2; % Events CAN be closer but we skip a lil ahead
            end
        end
        if ~isempty(event_vals)
            firedNeurons=[firedNeurons,i];
            firingThresholds=[firingThresholds,threshold];
            avgEventSizes=[avgEventSizes,sum(event_vals)/number_of_events];
            firingTimes = [firingTimes,cellFireTimes];
        end
    end
    %%% Prints out some info
    % DATA
    %  totalCells  activeCells  number_of_events  std_neg
    %  avgEventSizes()  event_vals()  intensityData(firedNeurons(i),:)  firedNeurons()
    %  std_negVals()  amplitudes() firingTimes()
    %  seconds  minutes
    
    %%% binaryFiring only keeps cells that have fired.
    binaryFiring(all(binaryFiring==0,2),:)=[];
    seconds = frames/fps;
    minutes = seconds/60;
    activeCells=length(firedNeurons);
    events_per_active_cell=number_of_events/activeCells;
    avgIntensity = mean(mean(intensityData));
    avgStdDev = mean(std_negVals);
    activePerMin = zeros(1,5); % holds percent active cells for each minute
    %%% This whole loop is for calculating % active cells per minute
    %%% You need to go minute by minute and average the result
    for index = 1:length(firingTimes)
        if firingTimes(index)>1
            subindex = index+1;
            singleActivePerMin = zeros(1,5);
            while firingTimes(subindex)<1 && subindex < length(firingTimes)
                frm = round(1/firingTimes(subindex));
                minuteIndex = ceil((frm / fps)/60); 
                singleActivePerMin(minuteIndex)=1; % will cause error if min>5
                subindex=subindex+1;
            end
            activePerMin = activePerMin + singleActivePerMin;
        end
    end
    %%% Now we combine the data and calculate percent active cells from the
    %%% active cells for each minute, last minute is weighted less than the
    %%% others
    avgActiveCellsPerMin = 0;
    for i=1:ceil(minutes)
        if i == ceil(minutes)
            % on last minute, average together the previous minutes and
            % then add the less-weighted last minute data
            avgActiveCellsPerMin = (avgActiveCellsPerMin+activePerMin(i)*(minutes-floor(minutes)))/minutes;
            percentActiveCells = 100 * avgActiveCellsPerMin / realCells;
        else
            % add together if not on last minute
            avgActiveCellsPerMin = avgActiveCellsPerMin + activePerMin(i);
        end
    end
 %   percentActiveWeighted = (events_per_active_cell*100*activeCells/totalCells)/minutes;
    
    % threshold_factor MINIMUM value is 7, otherwise breaks down
    % threshold_factor MAXIMUM value is 10, otherwise large data loss
    disp(' ')
    disp(outputName)
    disp(['  Threshold factor is ',num2str(threshold_factor)])
    disp(['  Movie is ',num2str(minutes),' minutes.'])
  %  disp(' ')
    % Mean intensity should be very close to zero (noise) but positive
    disp(['Mean intensity for all cells: ', num2str(avgIntensity)])
    disp([num2str(activeCells), ' active cells out of ', num2str(realCells), ' cells.'])
    disp(['Average negative std deviation: ', num2str(avgStdDev)])
    disp(['Lowest and Highest Neg.Std.Dev: ', num2str(min(std_negVals)), ', ', num2str(max(std_negVals))])
    disp(['Size of negative std deviation divided mean intensity: ',num2str(avgStdDev/avgIntensity)])
    disp(['Percent active cells: ',num2str((100*activeCells/realCells)),'%'])
    disp(['Percent active cells per minute: ',num2str(percentActiveCells),'%'])
    disp(['Number of events per minute: ',num2str((number_of_events)/minutes)])
    disp(['num events per minute per active cell: ',num2str((number_of_events/activeCells)/minutes)])
    disp(['Avg amplitude per event: ',num2str(mean(amplitudes))])
    disp(['Lowest and Highest amplitude: ', num2str(min(amplitudes)), ', ', num2str(max(amplitudes))])
    disp(['Avg area under curve per event: ',num2str(mean(avgEventSizes))])
    disp(' ')
    realCells
    numCellPlots = ceil(length(firedNeurons)/10);

    if length(firedNeurons)~=0
        Intensity_DeltaF_Plotter();
        saveas(gcf,outputName)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ ] = Intensity_DeltaF_Plotter()
    %INTENSITY_DELTAF_PLOTTER plots deltaF intensityData
    iStart = 1 + (cellPlot-1)*10;
    iEnd = min([iStart+9,length(firedNeurons)]);
    
    numPlots=length(firedNeurons(iStart:iEnd));

    % plotting the actual data
    hold off; % clear remaining figures
    clf
    inc = 1;
    for i=iStart:iEnd
        firedNeurons(i);
        plot( intensityData(firedNeurons(i),1:size(intensityData,2)-10)/1.25+(inc-1),'linewidth',0.5, 'color', 'k')
        inc = inc+1;        
        hold on;
        
    end
    % done plotting data

    % code for plotting a horozontal line
    %{m
    L_min = 0;
    L_max = length(intensityData(1,:));
    % L_min = L_max/12;
    y = 0; % constant y value
    d = 1; % direction
    r = 15;% y range
    a = linspace(L_min, L_max, 2);
    b = linspace(y, y, 2);
%     plot(a,b,'--k','linewidth',1.5)
    for i=1:((numPlots)*2-1)
        y = i/2;
        if mod(i,2)==1%%% If this horozontal line is a threshold
            y = round(y-0.5)+firingThresholds(round(i/2));
            a = linspace(L_min, L_max, 2);
            b = linspace(y, y, 2);
%             plot(a,b,'--k','linewidth',0)
        else%%% If line is an avg
            a = linspace(L_min, L_max, 2);
            b = linspace(y, y, 2);
%             plot(a,b,'--k','linewidth',0)
        end
    end
    %}
    % done plotting horozontal lines

    % Should eliminate white space from the plot (only works w/1 plot)
    %set(gca,'LooseInset',get(gca,'TightInset'))
    % Chooses range of the plot
    axis([0,L_max,-0.5,numPlots-0.2])
    % Sets the y tick labels
    set(gca,'YTick',(iStart-1):(iEnd-1),'YTickLabel',firedNeurons)
    set(gca,'YTick',0:numPlots-1,'YTickLabel',firedNeurons(iStart:iEnd) )
    set(gca, 'xtick', [-10000000 10000000]);
	set(gca, 'ytick', [-10000000 10000000]);
    title(['Cell Intensities with Firing Thresholds (factor=',num2str(threshold_factor),')'])
    ylabel('Cell Number')
    xlabel('Bottom: Frames, Top: Seconds')
    % This chunk will create a second set of Axis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %{
    ax1_pos = get(gca,'Position'); % position of first axes
    ax2 = axes('Position',ax1_pos,'XAxisLocation','top',...
        'YAxisLocation','right','Color','none');
    % sets the y tick labels
    set(ax2,'YTick',0:1,'YTickLabel',['' ''])
    seconds=round(length(intensityData(1,:))/3.91);
    set(ax2,'XTick',0:10,'XTickLabel',[0 seconds])
    %}
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function sliderCallback(src, evt)
  %  sliderVal=round(get(src, 'Value'));
    sliderVal=get(src, 'Value');
    fprintf('Slider value is: %d\n', sliderVal );
    threshold_factor = sliderVal;
    Calculate_Events(); 
    init_ui();
end

function init_ui()
    slmin=4.5;slmax=10;
    uicontrol('Style','slider','Callback',@sliderCallback,'Min',slmin,'Max',slmax,...
                    'SliderStep',[.5 .5]./(slmax-slmin),'Value',sliderVal,...
                    'Position',[5 5 200 20]);
end

pause;
  % DATA
    %  totalCells  activeCells  number_of_events  std_neg
    %  avgEventSizes()  event_vals()  intensityData(firedNeurons(i),:)  firedNeurons()
    %  std_negVals() firingThresholds()
    %  seconds  minutes
save(outputName,'intensityData', 'firedNeurons', 'number_of_events', 'threshold','amplitudes','avgEventSizes','firingThresholds','firingTimes','totalCells','binaryFiring');
close all

end