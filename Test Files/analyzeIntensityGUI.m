function [ ] = analyzeIntensityGUI( intavgNPCorr, outputName)
% analyzeIntensity converts a matrix into it's deltaF over avgF form,
% graphs the result, and analyzes the data
%% Defining initial variables
ignoreInitialSpiking = false; % Line 140 to change analyzable part
% Multiplied by negative standard deviation is the threshold for firing
% threshold_factor MINIMUM value is 7, otherwise breaks down
% threshold_factor MAXIMUM value is 10, otherwise large data loss
threshold_factor = 7.0; % threshold = std_neg * threshold_factor;
%%% Any especially noisy neurons above this threshold are labelled as not
%%% firing
stdCutoff = 0.053;
stdCutoff = 0.09;
stdCutoff = 1000;
stdCutoff = 0.099;

%%% Removes "cells" with an average intensity below a threshold
%%% indicating that it might be just dark neuropil
darknessThreshold = 600; % Above 1000 seems to strongly indicate a real cell
darknessThreshold = 0;
darknessThreshold = 0;


% Will run intensity data through a high pass filter
High_Pass_Filter = true;
Exponential_Fit = false;

initialTotalCells = length(intavgNPCorr(:,1));
std_negVals = ones(length(intavgNPCorr(:,1)),1);
fps = 3.91; % frames per second of the movie
removed_cells = zeros(length(intavgNPCorr(:,1)),1);

% This designates what kind of HP filter you use
hpFilt = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',0.00000000001,'PassbandRipple',0.1, ...
          'SampleRate',3.81);
%fvtool(hpFilt) % Displays the properties of the filter
%pause(2)

%% Indicates what cells are too dark to be considered cells
r=1;
while r < length(intavgNPCorr(:,1))
    darkLevel = mean(intavgNPCorr(r,:));
    if darkLevel < darknessThreshold
        removed_cells(r) = 1;
        fprintf('Removing a cell, too dark, light level: %d\n',darkLevel);
    end
    brightness_vals(r) = darkLevel;
    r=r+1;
end
fprintf('Avg cell brightness: %d\n', mean(brightness_vals));

%% Fitting function, Calculating Delta F over F, Indicates noisy cells
h = waitbar(0,'Finding DF/F and fitting to exponential for each cell');
%%% Intensity data is intavgNPCorr but with problematic cells removed
intensityData = intavgNPCorr;
frames = length(intavgNPCorr(1,:));
%%% Two term exponential fit to the data
x = linspace(1,frames,frames);
% Applies the exponential fit to the data and indicates noisy cells
r=1; % r===cell#
while r <= initialTotalCells
    if Exponential_Fit
        %%% Exponential fitting function
        cellIntensity = intensityData(r,:);
        x = (1:1:size(cellIntensity'));
        f = fit(x',cellIntensity','exp2','Start',[1000,-0.001,3500,0.0001]);
        for c=1:frames % Calculates DeltaF/F using exponential fit as avg
             intensityData(r,c) = (intensityData(r,c)-f(c))/f(c);
        end
    else
        %intensityData(r,:) = 2*intensityData(r,:)/max(intensityData(r,:));
    end
    %%% Forces the average to be 0
    intensityData(r,:)=intensityData(r,:)-mean(intensityData(r,:));
    %%% Remove all "cells" with particularly large standard deviations
    A = intensityData(r,:);% For some reason 'A' is needed
    std_neg = std(A(A<0));% Negative Std Dev of the signal
    std_negVals(r) = std_neg;
    if std_neg>stdCutoff
        removed_cells(r) = 1;
        fprintf('Removing a cell, too noisy, -std dev: %d\n',std_neg);
    end
    r=r+1;
    waitbar(r/initialTotalCells)
end
close(h)
fprintf('Avg brightness negative std deviation: %d\n', mean(std_negVals));

%% Applies a high pass filter onto each cell
%%% Test which filter is most optimal with /.TestFiles/FilterTestRealData.m
if High_Pass_Filter
    h = waitbar(0,'Putting each cell signal through a High Pass Filter');
    % hpfilt defined earlier, holds filter information
     fvtool(hpFilt) % Displays the properties of the filter
    for r=1:initialTotalCells % r===cell#
        %dataMir = [-fliplr(intensityData(r,:)), intensityData(r,:), -fliplr(intensityData(r,:))];
        %dataMir = filter(hpFilt,dataMir);
        %size(dataMir(floor(frames):floor(frames*2)-1))
        %size(intensityData(r,:))
        %intensityData(r,:) = dataMir(floor(frames):floor(frames*2)-1);
        intensityData(r,:) = filter(hpFilt,intensityData(r,:));
        waitbar(r/initialTotalCells)
    end
    close(h)
end

%% Removes cells marked as unusable or fake
%%% Intensity data is intavgNPCorr but with problematic cells removed
intensityData = intensityData(~removed_cells,:);
totalCells = length(intensityData(:,1));
std_negVals = std_negVals(~removed_cells);
brightness_vals = ones(length(intavgNPCorr(:,1)),1);
fprintf('Removed %d cells.\n', length(intavgNPCorr(:,1))-totalCells);

%% Defining variables used for calculations and initializes UI
assignin('base', 'Intensity_DeltaF', intensityData);
firedNeurons=[]; % An array that holds each firing cell's number
firingThresholds=[]; % The minimum threshold the signal must reach
firingTimes=[]; % Holds the frame of every cell that fires
%%% firingTimes is weird. It holds both the cell numbers and their firing times.
%%% #>1 is cell number, followed by (1/frameFired). All numbers greater
%%% than one is a cell number and every number less than one is the
%%% inverted frame on which it fired
% BinaryFiring is a matrix, cell activity on a frame is marked with "1"
binaryFiring = zeros(totalCells,frames);
cellPlot = 1;
numCellPlots = 1;

f = figure('KeyPressFcn',@keypress);
Calculate_Events();
slmin=4.5;slmax=10;
uicontrol('Style','slider','Callback',@sliderCallback,'Min',slmin,'Max',slmax,...
                    'SliderStep',[.5 .5]./(slmax-slmin),'Value',threshold_factor,...
                    'Position',[5 5 200 20]);
sliderVal=threshold_factor;

%% Handles every keypress
function keypress(src, evt)
    val = double(get(f,'CurrentCharacter'));
    switch val%%%% REPLACE THIS CRAP WITH ANOTHER SLIDER THAT CHOOSES WHICH PLOT
        case 28 % <-
            if cellPlot <= 1
                fprintf('Already on first plot\n');
                cellPlot=1;
            else
                cellPlot=cellPlot-1;
                Intensity_DeltaF_Plotter();
            end
            init_ui();
        case 29 % ->
            if cellPlot >= numCellPlots
                fprintf('Already on last plot\n');
            else
                cellPlot=cellPlot+1;
                Intensity_DeltaF_Plotter();
            end
            init_ui();
        otherwise
            fprintf('Button not recognized, value is: %d\n',val);
    end
    
end

%% This function iterates through every point looking for firing events 
function Calculate_Events()
    %%% Initializing variables to be used for each cell
    amplitudes = [];
    number_of_events = 0;
    avgEventSizes = [];
%%% For every segmented cell (cell i) run intensity analysis
%%% cell i
    for i=1:length(intensityData(:,1))%%% Data analysis
        cellFireTimes = (i); % Keeps track of the times at which the cell fires
        A = intensityData(i,:);
        threshold = std_negVals(i)*threshold_factor;
        event_vals=[]; % An array holding every single event size
%%%[FOR THE LAST CELL THIS CODE RUNS IMPROPERLY FOR SOME REASON]        
        startFrame = 5;
        if ignoreInitialSpiking
            startFrame = round(length(A)/8);
        end
        %%% Registers an event only if a value above the threshold is preceded by a value below it
        %%% frame ii
        for ii=startFrame:(length(A)-2) 
            %%% The following if statement will determine if a neuron has fired
            if (A(ii)<threshold) && (A(ii+1)>threshold) && (A(ii+2)>0) && std_negVals(i)<stdCutoff
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
         %       ii=ii+round(.8*fps);
                 ii=ii+2; % Events CAN be closer but we skip a lil ahead
            end
        end
        if ~isempty(event_vals)
            firedNeurons=[firedNeurons,i];
            firingThresholds=[firingThresholds,threshold];
            avgEventSizes=[avgEventSizes,sum(event_vals)/number_of_events];
            firingTimes = [firingTimes,cellFireTimes];
        end
        if isempty(event_vals)
            firingThresholds=[firingThresholds,0];
        end
    end
    %% Calculating real data 
    %%% Prints out some info
    %       currently "Active Cell" has fired at least once over course of video
    %       definition of active cell should be changeable
    % DATA
    % 0th Order: Total Cells, Active Cells, Total Events
    % 1st Order: Percent Active Cells, Normalized Events per Active/Total Cells
    %               Sum (# Events/Min) over Active Cells
    % 2nd Order: Synchrony, lower bounded events/(active*min)
    
    % [old]
    %  totalCells  activeCells  number_of_events  std_neg
    %  avgEventSizes()  event_vals()  intensityData(firedNeurons(i),:)  firedNeurons()
    %  std_negVals()  amplitudes() firingTimes()
    %  seconds  minutes
    
    % number_of_events
    
    %%% binaryFiringAll keeps inactive cells
    binaryFiringAll = binaryFiring;
    %%% binaryFiring only keeps cells that have fired.
    binaryFiring(all(binaryFiring==0,2),:)=[];
    seconds = frames/fps;
    minutes = seconds/60;
    activeCells=length(firedNeurons);
    events_per_active_cell=number_of_events/activeCells;
    
    totalEvents = length(firingTimes(firingTimes<1));
    avgStdDev = mean(std_negVals);
    avgIntensity = mean(mean(intensityData));

    % Variables filled inside the following loop
    activeCellsEachMin = zeros(1,ceil(minutes));
    eventsEachMin = zeros(1,ceil(minutes));
    %%% This loop is for calculating % active cells per minute
    index = 1; %%% You need to go minute by minute and average the result
    while index < length(firingTimes) % Goes through 'firingTimes'
        if firingTimes(index) >= 1 % if it contains a number >1 than it's the cell #
            cellNum = firingTimes(index);
            singleActivePerMin = zeros(1,ceil(minutes)); % Binary flag if cell was active during a minute
            subindex = index+1;
            while firingTimes(subindex)<1 && subindex<length(firingTimes)
                fireTime = round(1/firingTimes(subindex));
                minuteIndex = ceil((fireTime / fps)/60);
                singleActivePerMin(minuteIndex) = 1; % Marks cell being active during this minute
                eventsEachMin(minuteIndex) = eventsEachMin(minuteIndex) + 1; % Adds each event
                subindex=subindex+1;
            end
            activeCellsEachMin = activeCellsEachMin + singleActivePerMin;
            index = subindex;
        end
    end
    eventsEachMin;      %%% Calculated, binned by minute
    activeCellsEachMin; %%% Calculated, binned by minute
    
    leftoverMinute = minutes - floor(minutes); % what fraction of a minute is left over
    
    %%% Now we combine the data and calculate percent active cells from the
    %%% active cells for each minute, last minute is weighted less than the
    avgActiveCellsPerMin = 0; %%%  others based on leftoverMinute
    for i=1:ceil(minutes)
        if i == ceil(minutes)
            % on last minute, average together the previous minutes and
            % then add the less-weighted last minute data
            avgActiveCellsPerMin = (avgActiveCellsPerMin+activeCellsEachMin(i)*(leftoverMinute))/minutes;
            percentActiveCells = 100 * avgActiveCellsPerMin / totalCells;
        else
            % add together if not on last minute
            avgActiveCellsPerMin = avgActiveCellsPerMin + activeCellsEachMin(i);
        end
    end
    % avgActiveCellsPerMin bins minute by minute and so is always a lower bound
    avgActiveCellsPerMin;   %%% Calculated
    
    %%% 
    amplitudes
    avgEventSizes
    pause()
    
    
    %%% Printing out the calculated data
     disp(' ')
     disp(' ')
     disp(outputName)
     disp(['  Threshold factor is ',num2str(threshold_factor)])
     disp(['  Movie is ',num2str(minutes),' minutes.'])
    if seconds<180
        disp('Movies less than 180 seconds may not be representative of actual neural behavior')
    end
    % Mean intensity should be very close to zero (noise) but positive
    disp(' ')
     disp('ZEROTH ORDER DATA')
     disp(['Number of events total: ',num2str(number_of_events)])
     disp(['Number of total cells:  ',num2str(totalCells)])
     disp(['Number of Active cells: ',num2str(activeCells)])
     disp(['Percent Active Cells:   ',num2str(100*activeCells/totalCells),'%'])
     disp(' ')
     disp('FIRST ORDER DATA')
     disp(['Events per detected cell: ',num2str(number_of_events/totalCells)])
     disp(['Events per active cell:   ',num2str(number_of_events/activeCells)])
     disp(['Events per detected cell per minute: ',num2str(number_of_events/(totalCells*minutes))])
     disp(['Events per active cell per minute:   ',num2str(number_of_events/(activeCells*minutes))])
     disp(' ')
     disp(['Mean intensity for all cells: ', num2str(avgIntensity)])
     disp([num2str(activeCells), ' active cells out of ', num2str(totalCells), ' cells.'])
     disp(['Average negative std deviation: ', num2str(avgStdDev)])
     disp(['Lowest and Highest Neg.Std.Dev: ', num2str(min(std_negVals)), ', ', num2str(max(std_negVals))])
     disp(['Size of negative std deviation divided mean intensity: ',num2str(avgStdDev/avgIntensity)])
     disp(['Percent active cells: ',num2str((100*activeCells/totalCells)),'%'])
     disp(['Number of events per minute: ',num2str((number_of_events)/minutes)])
    %%% FIX THIS ONE
     disp(['num events per minute per active cell: ',num2str((number_of_events/activeCells)/minutes)])
     disp(['[NOT STATIONARITY] Percent active cells per minute (lower bound): ',num2str(percentActiveCells),'%'])
     disp(['[NOT STATIONARITY] Avg active cells per minute (lower bound): ',num2str(avgActiveCellsPerMin)])
     disp(['Avg amplitude per event: ',num2str(mean(amplitudes))])
     disp(['Lowest and Highest amplitude: ', num2str(min(amplitudes)), ', ', num2str(max(amplitudes))])
     disp(['Avg area under curve per event: ',num2str(mean(avgEventSizes))])
     disp(' ')
    
    numCellPlots = ceil(length(firedNeurons)/10);
    
    %{
    meandata = zeros(frames);
    for p=1:frames
        meandata(p) = mean(intensityData(:,p));
    end
    plot(meandata)
    pause()
    %}
    if length(firedNeurons)~=0
        Intensity_DeltaF_Plotter();
        saveas(gcf,outputName)
    end
end


%% Plotting code
function Intensity_DeltaF_Plotter()
    %INTENSITY_DELTAF_PLOTTER plots deltaF intensityData
    iStart = 1 + (cellPlot-1)*10;
    iEnd = min([iStart+9,length(firedNeurons)]);
    
    numPlots=length(firedNeurons(iStart:iEnd));

    % plotting the actual data
    hold off; % clear remaining figures
    clf
    inc = 1;
    for i=iStart:iEnd
        plot( intensityData(firedNeurons(i),:)+(inc-1),'linewidth',1.1)
        inc = inc+1;
        hold on;
    end
    % done plotting data

    % code for plotting a horozontal line
    L_min = 0;
    L_max = length(intensityData(1,:));
    % L_min = L_max/12;
    y = 0; % constant y value
    d = 1; % direction
    r = 15;% y range
    a = linspace(L_min, L_max, 2);
    b = linspace(y, y, 2);
    plot(a,b,'--k','linewidth',1.5)
    for i=1:((numPlots)*2-1)
        y = i/2;
        if mod(i,2)==1%%% If this horozontal line is a threshold
            %y = round(y-0.5)+firingThresholds(round(i/2));
            neuronNum = firedNeurons(iStart+round(i/2)-1);
            y = (i-1)/2+firingThresholds(neuronNum);
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
    axis([0,L_max,-0.5,numPlots-0.2])
    % Sets the y tick labels
%    set(gca,'YTick',(iStart-1):(iEnd-1),'YTickLabel',firedNeurons)
    set(gca,'YTick',0:numPlots-1,'YTickLabel',firedNeurons(iStart:iEnd) )
    title(['Cell Intensities with Firing Thresholds (factor=',num2str(threshold_factor),')'])
    ylabel('Cell Number')
    xlabel('Bottom: Frames, Top: Seconds')
    % This chunk will create a second set of Axis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ax1_pos = get(gca,'Position'); % position of first axes
    ax2 = axes('Position',ax1_pos,'XAxisLocation','top',...
        'YAxisLocation','right','Color','none');
    % sets the y tick labels
    set(ax2,'YTick',0:1,'YTickLabel',['' ''])
    seconds=round(length(intensityData(1,:))/3.91);
    set(ax2,'XTick',0:10,'XTickLabel',[0 seconds])
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
save(outputName,'intensityData', 'firedNeurons', 'number_of_events', 'threshold',...
    'amplitudes','avgEventSizes','firingThresholds','firingTimes','totalCells','binaryFiring');
save('stddata','std_negVals')
close all

end