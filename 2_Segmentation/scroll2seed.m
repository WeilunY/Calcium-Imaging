function [corrTHR, newSeeds]=scroll2seed(corrMat, image, neuronSize)
% 

% INSTRUCTIONS:
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
%       'up arrow' - Decrease threshold
%       'down -'   - Increase threshold

refCircle = false;

xManual=[];
yManual=[];
newSeeds=[];

width = length(image(1,:,1));
height = length(image(:,1,1));

THR=mean2(corrMat)+2*std2(corrMat);
step=1.05;
THR = THR/step;
indexy=find(corrMat<THR);
corrTHR=corrMat;
corrTHR(indexy)=0;
rgb(:,:,2)=image;
rgb(:,:,3)=image;
rgb(:,:,1)=image+uint8(255*corrTHR);
warning off all;
figure('Units', 'Normalize', 'Name', 'Image <-> template correlation matrix')
   imagesc(corrMat)
f = figure('WindowButtonDownFcn',@manualPick,'KeyPressFcn',@keypress,...
            'Name', 'Adjust specificity x sensitivity ratio using the up and down arrow keys');
h=imshow(uint8(rgb),'Border','tight'); 
     
figureAxis=get(h);
init_ui();

%%% Need to resize stack to be as big as image using ImageBiggerer.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
%
function keypress(src, evt)
    val = double(get(f,'CurrentCharacter'));
    switch val
        case 32 % ' '
        case 28 % <-
        case 29 % ->
        case 31 % down arrow
            THR=THR*step;
            indexy=find(corrMat<THR);
            corrTHR=corrMat;
            corrTHR(indexy)=0;
            rgb(:,:,1)=image+uint8(255*corrTHR);
            imshow(uint8(rgb), 'Parent', figureAxis.Parent);
        case 30 % up arrow
            THR=THR/step;
            indexy=find(corrMat<THR);
            corrTHR=corrMat;
            corrTHR(indexy)=0;
            rgb(:,:,1)=image+uint8(255*corrTHR);
            imshow(uint8(rgb), 'Parent', figureAxis.Parent);
        case 114 % r
        case 112 % p
        case 116 % t
            refCircle = ~refCircle;
        case 84 % PC t
            refCircle = ~refCircle;
        otherwise
            fprintf('Button not recognized, value is: %d\n',val);
            stop = true;
    end
    init_ui();
end
function init_ui()
    imagesc(corrMat)
    h=imshow(uint8(rgb),'Border','tight'); 
    figureAxis=get(h);
    
    rgb(:,:,1)=image+uint8(255*corrTHR);
    rgb(:,:,2)=image;
    rgb(:,:,3)=image;
    
    viscircles([xManual yManual], 1*ones(1,length(xManual))',...
                                    'EdgeColor', 'b','LineWidth',4);
%     fprintf('Currently on frame: %d\n',i);
    % Prints a reference circle that shows approximately what size neurons
    %   will be looked for
    if refCircle
        % Shows Avg neuron size     *3.8883 (/2 ?) pixels/micron
        hold on
        rsmall = neuronSize(1)*3.888/2; rbig = neuronSize(2)*3.888/2;
        r = (rsmall+rbig)/2;
        th = 0:pi/50:2.05*pi;
        xunit = r * cos(th) + 50;
        yunit = r * sin(th) + 50;
        plot(xunit, yunit, 'LineWidth',2*(rbig-rsmall)/3.888,'Color','w');
        
        txt1 = 'Average Circle Size';
        text( 30, 10, txt1, 'color', 'w')
        
        hold off
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%

function manualPick(~,~)
    buttonTMP=1;
    while buttonTMP~=2
        imshow(uint8(rgb), 'Parent', figureAxis.Parent);
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
        
    end
    newSeeds=[xManual yManual];
%     init_ui();
%corrTHR(sub2ind(size(corrTHR), round(yManual), round(xManual)))=1;
end

pause();

close all %%%Comment out this "close all" to see two images - 
          %%%the one with neurons and correlation seeds and the image with 
          %%%neurons with overlying donuts.
end