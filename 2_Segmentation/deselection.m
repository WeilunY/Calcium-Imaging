function [maskStructure0,rgb0,maskStack0]=deselection(maskStructure, maskStack, rgb)
% Click and deselect cells
rgb0=rgb;
maskStructure0=maskStructure;
maskStack0=maskStack;
blueness = 255; % Determines how bright the selected cells are displayed

warning off all;
f = figure('WindowButtonDownFcn',@manualPick,'KeyPressFcn',@keypress,...
            'Name', 'Click anywhere, the nearest cell will be removed, press "r" to reset');
h=imshow(uint8(rgb0),'Border','tight'); 
% figureAxis=get(h);

% Initializes the blue circles to their values as originally calculated
init_ui();

%%% Need to resize stack to be as big as image using ImageBiggerer.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function keypress(src, evt)
    val = double(get(f,'CurrentCharacter'));
    switch val
        case 32 % ' '
            init_ui();
        %%% Reset button
        case 114 % r
            rgb0=rgb;
            maskStructure0=maskStructure;
            maskStack0=maskStack;
            fprintf('Wiped changes made to cell selection\n');
            init_ui();
        case 30 % up arrow
            blueness = blueness+25;
            if blueness>255
                blueness=255;
            end
            rgb0=rgb;
            init_ui();
            pause(0.01);
        case 31 % down arrow
            blueness = blueness-25;
            if blueness<5
                blueness=5;
            end
            rgb0=rgb;
            init_ui();
            pause(0.01);
        otherwise
            fprintf('Button not recognized, value is: %d\n',val);
    end
    
end
function init_ui()
    imS=size(maskStack0,1);
    %%% This generates the blue circles, the information of the circles is
    %%% held in maskStructure & maskStructure0
    for sindex=1:size(maskStructure0,2)
        xCent = maskStructure0(sindex).CenterCoor(1);
        yCent = maskStructure0(sindex).CenterCoor(2);
        rgb0(xCent-imS/2:xCent+imS/2-1, yCent-imS/2:yCent+imS/2-1, 3)=...
        rgb0(xCent-imS/2:xCent+imS/2-1, yCent-imS/2:yCent+imS/2-1, 3)+blueness*maskStack0(:,:,sindex);
    end
    h=imshow(uint8(rgb0),'Border','tight');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function manualPick(~,~)
    buttonTMP=1;
    while buttonTMP~=2
        [xTMP,yTMP,buttonTMP] = ginput(1);
        minDist = 9999;
        deleteIndex = -1;
        % Find the circle closest to where you click
        for i=1:size(maskStructure0,2)
               xDist = xTMP - maskStructure0(i).CenterCoor(2);
               yDist = yTMP - maskStructure0(i).CenterCoor(1);
               currDist = sqrt(xDist*xDist+yDist*yDist);
               if currDist < minDist
                   minDist = currDist;
                   deleteIndex = i;
               end
        end
        
        % Remove the closest circle by deleting the index of all variable clones
        if buttonTMP==1
            rgb0=rgb;
            maskStructure0(deleteIndex)=[];
            maskStack0(:,:,deleteIndex)=[];
            init_ui();
        elseif buttonTMP==3
        end
    end
%    init_ui();
end

pause();
close all %%%Comment out this "close all" to see two images - 
          %%%the one with neurons and correlation seeds and the image with 
          %%%neurons with overlying donuts.
end