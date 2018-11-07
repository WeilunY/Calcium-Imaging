function [maxHorizontalMovement, adjustedImage]=rowShifter2(rawImage,reference,maxHorizontalMovement, rByRcorr, penaltyPow)
reference=double(reference);
rawImage=double(rawImage);
penaltyParabolCoef=-1*rByRcorr*penaltyPow;%%%%Depends on average rowByRow correlation
%maxHorizontalMovement=10; % maximal displacement (from the reference) movements. In pixels;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

corrCoefBuffer=zeros(1,2*maxHorizontalMovement+1);
%displacementN=length(corrCoefBuffer);
displacement=-maxHorizontalMovement:maxHorizontalMovement;
adjustedImage=ones(size(rawImage), 'uint16');
adjustedImage=uint16(adjustedImage*mean2(reference));
penalty=polyval([penaltyParabolCoef 0 0], -2*maxHorizontalMovement:2*maxHorizontalMovement);
penaltyVertex=round(numel(penalty)/2);
maxMoveOLD=maxHorizontalMovement+1;

for index=maxHorizontalMovement+1:size(rawImage,1)-maxHorizontalMovement

    for dindex=-maxHorizontalMovement:maxHorizontalMovement %displacement index
        %ccTMP=corrcoef(double(reference(index,maxHorizontalMovement+1:end-
        %maxHorizontalMovement)), double(rawImage(index,maxHorizontalMovement+dindex+1:end-maxHorizontalMovement+dindex)));
        d_a=reference(index,maxHorizontalMovement+1:end-maxHorizontalMovement);
        d_b=rawImage(index,maxHorizontalMovement+dindex+1:end-maxHorizontalMovement+dindex);
        corrCoefBuffer(dindex+maxHorizontalMovement+1)=sum(sum(d_a.*d_b))/sqrt(sum(sum(d_a.*d_a))*sum(sum(d_b.*d_b)));
    end

%     maxMoveOLD
%     penalty
%     corrCoefBuffer
%     
% pause;
if mean(isnan(corrCoefBuffer))==1
   corrCoefBuffer=zeros(size(corrCoefBuffer));
end
if mean(isnan(corrCoefBuffer))>0
   %corrCoefBuffer
   corrCoefBuffer=zeros(size(corrCoefBuffer));
end
    corrCoefBuffer=corrCoefBuffer+penalty(penaltyVertex+1-maxMoveOLD:penaltyVertex+1+2*maxHorizontalMovement-maxMoveOLD);
    
    [~,maxMove]=max(corrCoefBuffer);%%%possibly to use a parabolical fit???
    maxMoveOLD=maxMove;
    %maxMove=displacementN+1-maxMove;
    adjustedImage(index,maxHorizontalMovement+1-displacement(maxMove):...
        end-maxHorizontalMovement-displacement(maxMove))=...
        uint16(rawImage(index,maxHorizontalMovement+1:end-maxHorizontalMovement));
    %  pause

end

% subplot(3,1,1)
% imshow(rawImage)
% subplot(3,1,2)
% imshow(adjustedImage)
% subplot(3,1,3)
% imshow(reference)
% mean2(abs(reference(:,2*maxHorizontalMovement+1:end-2*maxHorizontalMovement)
  %-rawImage(:,2*maxHorizontalMovement+1:end-2*maxHorizontalMovement)))
% mean2(abs(reference(:,2*maxHorizontalMovement+1:end-2*maxHorizontalMovement)
  %-adjustedImage(:,2*maxHorizontalMovement+1:end-2*maxHorizontalMovement)))
% pause;

end
