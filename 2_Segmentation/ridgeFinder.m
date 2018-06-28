function [ridgeRows]=ridgeFinder(polarImage, springCoeff)

%drawing=polarImage;
radiusMean=mean(polarImage,2);
[normMax,memberMax] = max(radiusMean);
radiusMeanN=radiusMean/normMax;
memberDiff2=diff(radiusMeanN,2);
%plot(memberDiff2)
[~, lowestDiff2Position] = min(memberDiff2);
springEq=round((memberMax+lowestDiff2Position+1)/2); %%%%%%Here the equilibrium 
%value is constructed from maximal intensity position in averaged-over-angle 
%image and from minimum second derivation from this curve
[rowN, angleNums]=size(polarImage);

ridgeRows=zeros(1,angleNums);
disconnected=1;
index=0;
rowTMP=springEq;
% runN=1;
while disconnected
index=mod(index,angleNums)+1;
if (rowTMP>1) && (rowTMP<rowN)
[~,stepCoeff] = max(polarImage(rowTMP-1:rowTMP+1,index));
rowTMP=(rowTMP+stepCoeff-2)+round(springCoeff*(springEq-(rowTMP+stepCoeff-2))); %%%The second term serves as a spring to stabilize the algorithm. Could be changed to higher order function in distance-equilibrium value.
elseif rowTMP==1
[~,stepCoeff] = max(polarImage(rowTMP:rowTMP+1,index));  
rowTMP=(rowTMP+stepCoeff-1)+round(springCoeff*(springEq-(rowTMP+stepCoeff-1))); %%%The second term serves as a spring to stabilize the algorithm. Could be changed to higher order function in distance-equilibrium value.
elseif rowTMP==rowN
[~,stepCoeff] = max(polarImage(rowTMP-1:rowTMP,index));  
rowTMP=(rowTMP+stepCoeff-2)+round(springCoeff*(springEq-(rowTMP+stepCoeff-2))); %%%The second term serves as a spring to stabilize the algorithm. Could be changed to higher order function in distance-equilibrium value.
end

disconnected=ridgeRows(index)~=rowTMP;
ridgeRows(index)=rowTMP;

%drawing(rowTMP,index)=255;

%pause;
% runN=runN+1;

end

% runN
% subplot(2,2,1);
%  imshow(drawing,'InitialMagnification',1000);
%  pause;
end