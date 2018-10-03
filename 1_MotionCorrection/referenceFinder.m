function [reference]=referenceFinder(stack, fps, startTime, calmTimeSt, calmTimeEnd, Nimages, percentile)
imageN=size(stack,3);
zeroPosition=ceil((calmTimeSt-startTime)*fps);
realNimages=floor((calmTimeEnd-calmTimeSt)*fps);

sureN=min(imageN-zeroPosition,realNimages); %%%This line to cope with an
%extremely unprobably, but yet possible situation, in which calm period is
%at the end of the record and fps number is not preciselly defined. For
%example stack "20160218_mouse5_001.tif" can cause such fault
step=floor(sureN/Nimages);
if step<1 && imageN-zeroPosition>0
    disp('Too short calm period. Increase the imdilate parameter in globalMC');
    return;
end
%%%%%%%%%%

corrMatix=zeros(Nimages);
h = waitbar(0,'Calculating reference image.');
size(stack)
Nimages
for index=1:Nimages
    waitbar(index/ Nimages)
    for index2=(index+1):Nimages
        index2;
        corrMatix(index, index2)=corr2(stack(:,:,zeroPosition+step*index), stack(:,:,zeroPosition+step*index2));
    end
end
close(h)

maxima=max(corrMatix);
%THR=prctile(maxima, 95);
sortedMaxima=sort(maxima); %%%% This line and the following one are faster than "prctile" function from Statistics Toolbox from Mathworks
THR=sortedMaxima(round(length(sortedMaxima)*percentile));

slices=find(maxima>THR);
reference=mean(stack(:,:,zeroPosition+step*slices),3);
reference=uint16(reference);

size(reference)

end