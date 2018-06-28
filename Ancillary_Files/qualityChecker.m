function [newSeeds,maskStructure0,rgb0,maskStack0]=qualityChecker( imageRef, stack, maskStructure, maskStack)
% Original Reference Image saved, 'bigStack' is the stack of images that
% has been enlarged to be viewed more easily

% image: whatever is currently the background
% imageRef: the initial
% var  : Initial state of inputted parameters, not changed
% var0 : Dynamic state of inputted parameters, user changes these

% maskStructure: struct array, each index has info on a different cell


global stop
global looping
stop = true;
looping = true;
deselectionMode = false;
instructions = true;
%%% image is the backgruond image
%%% rgb has the blue circles on it already

image = imageRef;
rgb = cat(3, image, image, image);
rgb0 = cat(3, image, image, image);
maskStructure0=maskStructure;
maskStack0=maskStack;

width = length(image(1,:,1));
height = length(image(:,1,1));
imageRef = image;
stack = stack(:,:,1:5);
stackLen = length(stack(1,1,:));
bigStack = uint8(zeros(512,2333,stackLen));
i = 1;
while i <= stackLen
    bigStack(:,:,i) = uint8(ImageBiggerer(stack(:,:,i)));
    fprintf('Stack enlarging, Image: %d out of %d\n', i, stackLen);
    i = i + 1;
end
i=1;
% These three variables record seeds the user wants to add
xManual=[];
yManual=[];
newSeeds=[];

step=1.5; % How much the color multiplier changes with each up/down arrow press
colorMult=200; % Out of 255, dictates how bright the red/blue circles are

warning off all;
f = figure('WindowButtonDownFcn',@manualPick,'KeyPressFcn',@keypress,...
            'Name', 'manual selection and deselection');
h=imshow(uint8(rgb0),'Border','tight');
figureAxis=get(h);

% Initializes the blue circles to their values as originally calculated
init_ui();

%%% Need to resize stack to be as big as image using ImageBiggerer.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
function keypress(src, evt)
    val = double(get(f,'CurrentCharacter'));
    switch val
        case 32 % ' '
            stop = ~stop;
            while i < length(bigStack(1,1,:)) && ~stop
                image = bigStack(:,:,i);
                init_ui();
                pause(0.15);
                if stop
                    break;
                end
                i = i+1;
                if i == stackLen && looping
                    i = 1;
                end
            end
            fprintf('Currently on frame: %d\n',i);
            init_ui();
        case 28 % <-
            stop = true;
            i = i-1;
            if i<1
                i = 1;
            end
            image = bigStack(:,:,i);
            fprintf('Currently on frame: %d\n',i);
            init_ui();
            pause(0.25);
        case 29 % ->
            stop = true;
            i = i+1;
            if i>stackLen
                i = stackLen;
            end
            image = bigStack(:,:,i);
            fprintf('Currently on frame: %d\n',i);
            init_ui();
            pause(0.25);
        case 30 % up arrow
            colorMult = colorMult*step;
            if colorMult > 255
                colorMult = 255;
            end
            init_ui();
            pause(0.2);
        case 31 % down arrow
            colorMult = colorMult/step;
            if colorMult < 10
                colorMult = 10;
            end
            init_ui();
            pause(0.2);
        case 105 % i [toggle instructions]
            instructions = ~instructions;
            init_ui();
            pause(0.2);
        case 114 % r [reference image]
            image = imageRef;
            i = 1;
            fprintf('Currently on Reference Image\n');
            init_ui();
        case 112 % p [pause]
            stop = true;
        case 116 % t [toggle mode]
            deselectionMode = ~deselectionMode;
            stop = true;
            init_ui();
        case 117 % u [undo all]
            maskStructure0=maskStructure;
            maskStack0=maskStack;
            newSeeds = [];
            xManual=[];
            yManual=[];
            image = imageRef;
            init_ui();
        otherwise
            fprintf('Button not recognized, value is: %d\n',val);
            stop = true;
    end
    
end
function init_ui()    
    if deselectionMode
       color = 1; % Circles in green if in deselection mode
    else
       color = 3; % Circles are blue if not in deselection mode
    end
    rgb0 = cat(3, image, image, image);

    %%% Following code is from deselection.m
    imS=size(maskStack0,1);
    %%% This generates the blue/red circles, the information of the circles is
    %%% held in maskStructure & maskStructure0
    for sindex=1:size(maskStructure0,2)
        xCent = maskStructure0(sindex).CenterCoor(1);
        yCent = maskStructure0(sindex).CenterCoor(2);
        rgb0(xCent-imS/2:xCent+imS/2-1, yCent-imS/2:yCent+imS/2-1, color)=...
         uint8(rgb0(xCent-imS/2:xCent+imS/2-1, yCent-imS/2:yCent+imS/2-1, color))...
          +uint8(colorMult*maskStack0(:,:,sindex));
    end
    h=imshow(uint8(rgb0),'Border','tight');
    displayIntructions()
    %%% End of deselection code
    
    % This displays all of the user-inputted cells
    viscircles([xManual yManual], 1*ones(1,length(xManual))',...
                                    'EdgeColor', 'b','LineWidth',4);
    
    %fprintf('Currently on frame: %d\n',i);
end
function displayIntructions()
    if instructions
    txt0 = 'Click anywhere to add/remove a cell (processed in next step)';
    txt1 = 'Arrow keys up/down change cell brightness';
    txt2 = 'Space to play/pause movie, arrow right/left selects frame';
    txt3 = 'p: emergegncy pause';
    txt4 = 'r: returns to reference image';
    txt5 = 't: toggles from blue selection mode to red deletion mode';
    txt6 = 'u: undo all changes and return to initial state';
    txt7 = 'i: toggle instructions';
    text( 10, 10, txt0, 'color', 'w')
    text( 10, 25, txt1, 'color', 'w')
    text( 10, 40, txt2, 'color', 'w')
    text( 10, 55, txt3, 'color', 'w')
    text( 10, 70, txt4, 'color', 'w')
    text( 10, 85, txt5, 'color', 'w')
    text( 10, 100, txt6, 'color', 'w')
    text( 10, 115, txt7, 'color', 'w')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%

%{
function figScroll(src,evnt)
    if evnt.VerticalScrollCount < 0 % scroll up
        colorMult = colorMult*step;
        if colorMult > 255
            colorMult = 255;
        end
    elseif evnt.VerticalScrollCount > 0 % scroll down 
        colorMult = colorMult/step;
        if colorMult < 10
            colorMult = 10;
        end
    end
    
    init_ui();
    pause(0.01);
end
%}

function manualPick(~,~)
    buttonTMP=1;
    while buttonTMP~=2
        % Clicking has a different function depending on what mode you're in 
%%% Selection Mode
        if ~deselectionMode % If you are not in deselection mode
            viscircles([xManual yManual], 1*ones(1,length(xManual))',...
                                        'EdgeColor', 'b','LineWidth',5);
            [xTMP,yTMP,buttonTMP] = ginput(1);
            if buttonTMP==1
                 xManual=[xManual;xTMP];
                 yManual=[yManual;yTMP];
                 %  DO NOT ADD if too close to the edge.
                 if xTMP<35 || xTMP>(width-35) || yTMP<35 || yTMP>(height-35)
                     xManual=xManual( 1 : length(xManual)-1 );
                     yManual=yManual( 1 : length(yManual)-1 );
                 end
            elseif buttonTMP==3
                rectTMP = getrect(f);    
                points2remove=find(xManual>rectTMP(1)&xManual<rectTMP(1)+...
                    rectTMP(3)&yManual>rectTMP(2)&yManual<rectTMP(2)+rectTMP(4));
                xManual(points2remove)=[];
                yManual(points2remove)=[];
            end
            newSeeds=[xManual yManual];
%%% Deselection Mode
        else% If you are in deselection mode
            [xTMP,yTMP,buttonTMP] = ginput(1);
            minDist = 9999;
            deleteIndex = -1;
            % Find the circle closest to where you click
            for j=1:size(maskStructure0,2)
                   xDist = xTMP - maskStructure0(j).CenterCoor(2);
                   yDist = yTMP - maskStructure0(j).CenterCoor(1);
                   currDist = sqrt(xDist*xDist+yDist*yDist);
                   if currDist < minDist
                       minDist = currDist;
                       deleteIndex = j;
                   end
            end
            % Remove the closest circle by deleting the index of all variable clones
            if buttonTMP==1
                size(maskStructure0);
                maskStructure0(deleteIndex)=[];
                maskStack0(:,:,deleteIndex)=[];
                size(maskStructure0);
                init_ui();
            elseif buttonTMP==3
            end
        end
        
    end
    init_ui();
%corrTHR(sub2ind(size(corrTHR), round(yManual), round(xManual)))=1;
end

pause;

close all %%%Comment out this "close all" to see two images - 
          %%%the one with neurons and correlation seeds and the image with 
          %%%neurons with overlying donuts.
end