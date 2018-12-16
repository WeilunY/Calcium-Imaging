function [corrTHR, newSeeds]=slider_plot(corrMat, image, stack)

global stop
global looping
stop = true;
looping = true;

width = length(image(1,:,1))
height = length(image(:,1,1))

imageRef = image;
stack = stack(:,:,1:20);
stackLen = length(stack(1,1,:));
bigStack = uint8(zeros(512,2333,stackLen));
i = 1;
while i <= stackLen
    bigStack(:,:,i) = uint8(ImageBiggerer(stack(:,:,i)));
    fprintf('Stack enlarging, Image: %d out of %d\n', i, stackLen);
    i = i + 1;
end
i=2;

xManual=[];
yManual=[];
newSeeds=[];

THR=mean2(corrMat)+2*std2(corrMat);
step=1.05;
indexy=find(corrMat<THR);
corrTHR=corrMat;
corrTHR(indexy)=0;
rgb(:,:,2)=image;
rgb(:,:,3)=image;
rgb(:,:,1)=image+uint8(255*corrTHR);
warning off all;
figure('Units', 'Normalize', 'Name', 'Image <-> template correlation matrix')
   imagesc(corrMat)
f = figure('WindowScrollWheelFcn',@figScroll, 'WindowButtonDownFcn',@manualPick,...
            'KeyPressFcn',@keypress,...
            'Name', 'Adjust specificity x sensitivity ratio using the mouse wheel');
h=imshow(uint8(rgb),'Border','tight'); 
     
figureAxis=get(h);

%%% Need to resize stack to be as big as image using ImageBiggerer.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%
%
function keypress(src, evt)
    val = double(get(f,'CurrentCharacter'));
    switch val
        case 32 % ' '
            stop = ~stop;
            while i < length(bigStack(1,1,:)) && ~stop
                image = bigStack(:,:,i);
                init_ui();
                pause(0.1);
                if stop
                    break;
                end
                i = i+1;
                if i == stackLen && looping
                    i = 1;
                end
            end
            init_ui();
        case 28 % <-
            stop = true;
            if i<=1
                i = 1;
            else
                i = i-1;
            end
            image = bigStack(:,:,i);
            init_ui();
            pause(0.25);
        case 29 % ->
            stop = true;
            if i>=stackLen-1
                i = stackLen-1;
            else
                i = i+1;
            end
            image = bigStack(:,:,i);
            init_ui();
            pause(0.25);
        case 114 % r
            image = imageRef;
            i = 1;
            fprintf('Currently on Reference Image\n');
            init_ui();
        case 112 % p
            stop = true;
        otherwise
            fprintf('Button not recognized, value is: %d\n',val);
            stop = true;
    end
    
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
    fprintf('Currently on frame: %d\n',i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%


function figScroll(src,evnt)
    if evnt.VerticalScrollCount > 0 
        %disp('ahoj1')
        THR=THR*step;
        indexy=find(corrMat<THR);
        corrTHR=corrMat;
        corrTHR(indexy)=0;
        rgb(:,:,1)=image+uint8(255*corrTHR);
        imshow(uint8(rgb), 'Parent', figureAxis.Parent);
    elseif evnt.VerticalScrollCount < 0 
        % disp('ahoj2')
        THR=THR/step;
        indexy=find(corrMat<THR);
        corrTHR=corrMat;
        corrTHR(indexy)=0;
        rgb(:,:,1)=image+uint8(255*corrTHR);
        imshow(uint8(rgb), 'Parent', figureAxis.Parent);
    end
 %   init_ui();
end

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
    init_ui();
%corrTHR(sub2ind(size(corrTHR), round(yManual), round(xManual)))=1;
end

pause;

close all %%%Comment out this "close all" to see two images - 
          %%%the one with neurons and correlation seeds and the image with 
          %%%neurons with overlying donuts.
end