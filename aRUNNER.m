function [] = aRUNNER()
% Runs Motion Correction and Segmentation algorithms together
%{
    %%%%% User Guide %%%%%
aRUNNER drives the entire process, it will run Motion Correction, then 
Segmentation, then finally some Intensity Analyses. This User Guide is
intended to make using this program easy, if any part is unclear or if you
ever encounter an error email me at adnewber@ucsd.edu.

Every file either needs to be in same the Matlab folder as this code (left)
or you need to specify the file's path.
Only a TIFF file is needed however a converted WCP files will make t1he
program run much, much better.
You will need:
    name.tiff
    name[1-1].mat [converted from name[1-1].wcp on the imaging computer]

MOTION CORRECTION
Once all the files are in this folder just run aRUNNER, Motion Correction
must run first and it will take at least 5 minutes to complete all the
loading bars. 
Once finished, the Motion Correlation Threshold figure will come up, the
user must choose a threshold that will cut off all parts of the movie where
the thick black line dips below the dotted threshold line. [needs work]

SEGMENTATION
first GUI
The first thing that comes up is a black and white avg image of the movie,
cells should be clearly visible. Use central mouse button to scroll up and
down, you'll see spots of red spread or recede. The red indicates that the
program thinks there is a cell in that area. Any area with clear cells and
no red you can click with the left mouse button, a small point will appear,
the program will add cells you select by hand as well as the ones it finds
automatically. Once finished type a character into the command window.
second GUI
Next is the deselection screen. The program will display in blue every cell
it found from you using the first GUI. Click anywhere and the nearest cell
will be removed. Use this to remove false positives.
third GUI
This final GUI is a combination of the first and second. It allows you to
select cells and remove existing ones, switch from selection mode to
removal mode by pressing 't'. In selection mode all the cells are blue and
it functions identical to the first GUI. In deletion mode it function
identical to the second GUI. Use the scrollwheel to change intensity of the
blue/red to see the background easier.

INTENSITY ANALYSIS
A plot of every cell's intensity over time will be displayed. Use left and
right arrow keys to see more. Doesn't need user intervention, time into
command window to move on.
More plots will come up displaying various information.


% Selection GUI INSTRUCTIONS:
%   Initially there will be a black and white reference image with
%   splotches of red which indicate areas found to resemble a donut
%   cell shape. SPACEBAR will play and pause the movie (stack) and you will
%   be able to manually select cells by clicking on the image/movie, this
%   brings you into cell selection mode. After you finish selecting cells 
%   you must press ENTER to leave cell selection mode. Press 'p' to just
%   pause the movie, may be easier to pause with 'p' than with spacebar at
%   times. Press 'r' to bring up the reference image again (the cleanest
%   image by far).
%       'SPACEBAR' - Play/Pause toggle
%       'ENTER'    - Exit cell selection
%       'p'        - Pause
%       'r'        - Reference Image
%       '<-'       - One frame backwards
%       '->'       - One frame forwards
%       't'        - Alternate between cell selection/removal [GUI 3]
    %%%%% User Guide %%%%%
%}

addpath('./1_MotionCorrection/');
addpath('./2_Segmentation/');
addpath('./3_DataExtraction/');
addpath('./4_DataAnalysis/');

% Run(MC/Seg/IntAnl) determines which parts of the program will run
RunMC  = 1;
RunSeg = 1;
RunIntAnl = 1;
quickIntensityAnalysis = 0; % Will set other options to 0

motionData = true; %%% Set to false if no WCP file
% if false then the motion correctioin will be less accurate unlesse
% you change calm period by hand
performMotionCorrection = true; %%% True by default, setting to false 
%%% prevents the code from running motion correction

region = 'DG'; % Which region of the brain is this movie in? (DG, CA1)
               % used in segmCorr to determine cell sizes
               % Pixels per micron: 3.8883 (/2 ?)
% region = | Hilus | DG | CA1 | CA3 |

%%% GUI allows parameter selection %%%
[motionData, region, RunMC, RunSeg, RunIntAnl, quickIntensityAnalysis,...
    FileName, PathName] = initialGUI();

%%% Creates a folder called 'data' in the directory we work with
if ~exist([PathName,'Data/'], 'dir') %% Check if there is a folder called 'data'
    %If folder does not exist, create it.
    mkdir(PathName, 'Data')
end

% This needs to match the name of your tiff file
file = FileName;
% Path to the tiff file
FileName1 = char(strcat(FileName,".tif")); 
% Path to the matlab file
FileName2 = char(strcat(FileName,"[1-1].mat")); % Convention: file.tiff => file[1-1].mat
file = char(strcat(PathName, file));
fileData = char(strcat(PathName, 'Data/', FileName));
PathName = char(PathName);

file %%% Full PATHNAME + FILENAME for the original files
fileData; %%% Full PATHNAME/Data + FILENAME for outputing data

%qualityChecker(imageCN, stackAdjustedCut, maskStructInt, maskStack, intensityData)

%%%%% Run Motion Correction
if RunMC && ~quickIntensityAnalysis 
    globalMC(FileName1,PathName,FileName2,PathName,motionData);
    load([fileData,'.mat']);
    if motionData
        load([PathName,FileName2]);
    else
        Y1 = 0;
    end
        %%%%% Run Motion Removal To Clean Up Movie
    motionCutterGUI(motionCompensation, stackAdjusted,Y1,[fileData,' motionCut_Output'],motionData);
    load([fileData,' motionCut_Output.mat']);
        %%%%% Save The Final Tiff Stack For Reference Later
    toTiff(stackAdjustedCut,[fileData,' CUT']);
end
    %%%%% Run Segmentation 
if RunSeg  && ~quickIntensityAnalysis 
    load([fileData,' motionCut_Output.mat']);
    tif2IntensitiesMatt(stackAdjustedCut,10,[fileData,' segOutput'], region);
end
    %%%%% Run Intensity Analysis
if quickIntensityAnalysis || RunIntAnl
    load([fileData,' motionCut_Output.mat']);
    load([fileData,' segOutput.mat']);
    analyzeIntensityGUI( intavgNPCorr, [fileData,' intensity data']);
    load([fileData,' intensity data.mat']);
    intensityZ( stackAdjustedCut, intavgNPCorr, firingTimes, binaryFiring, firedNeurons, fileData);
end


end