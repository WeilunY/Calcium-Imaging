function [maxXMovement, adjustedImage, M]=globalShifter2(rawImage, reference,maxDisplacement)
% Change data types to type double
reference=double(reference);
rawImage=double(rawImage);
%%% X and Y displacement might be independent if desirable.
maxXMovement=maxDisplacement;%10; % maximal X displacement (from the reference) movements. In pixels;
maxYMovement=maxDisplacement;%10; % maximal Y displacement (from the reference) movements. In pixels;

corrCoefBuffer=zeros(2*maxYMovement+1,2*maxXMovement+1);

displacement=-maxXMovement:maxXMovement;
adjustedImage=zeros(size(rawImage), 'uint16');

for yindex=-maxYMovement:maxYMovement %y displacement index
    for xindex=-maxXMovement:maxXMovement %x displacement index
        d_a=reference(maxYMovement+1:end-maxYMovement,maxXMovement+1:end-maxXMovement);
        d_b=rawImage(maxYMovement+yindex+1:end-maxYMovement+yindex,maxXMovement+xindex+1:end-maxXMovement+xindex);
        corrCoefBuffer(yindex+maxYMovement+1,xindex+maxXMovement+1)= ...
            sum(sum(d_a.*d_b))/sqrt(sum(sum(d_a.*d_a))*sum(sum(d_b.*d_b)));
    end
end
linearBuffer=corrCoefBuffer(:);
[M,I] = max(linearBuffer);
[I_row, I_col] = ind2sub([size(corrCoefBuffer,1) size(corrCoefBuffer,2)],I);

adjustedImage(maxXMovement+1-displacement(I_row):end-maxXMovement-displacement(I_row),...
    maxXMovement+1-displacement(I_col):end-maxXMovement-displacement(I_col))=uint16(...
    rawImage(maxYMovement+1:end-maxYMovement,maxXMovement+1:end-maxXMovement));


end
