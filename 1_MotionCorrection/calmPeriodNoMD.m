function [ calmTimeSt,calmTimeEnd, M_calm, M_mov ] = calmPeriodNoMD( stack )
% Calm Period finding when you do not have Motion Data.
frames = size(stack,3);

% Create an average over each frame, save as stackMean
stackSum = sum(stack,3);
stackMean = uint16(stackSum/size(stack,3));
%imshow(stackMean)

% Creates an array that holds correlation coefficients between stackMean
% and each frame of the stack
corrWithMean = zeros(frames,1);
for i = 1:frames
    corrWithMean(i) = corr2(stack(:,:,i),stackMean);
end

% Looking for a new, more stable reference image for more accurate
% correlation data. Slide a window across and find maximum chunk of
% correlation.
windowSize = 100;
maxWindowCorr = 0;
maxWindowCorrIndex = -1;

for i = 1:(frames-windowSize)
    currWindowCorr = sum( corrWithMean(i:i+10) );
    
    if currWindowCorr > maxWindowCorr
        maxWindowCorr = currWindowCorr;
        maxWindowCorrIndex = i;
    end
end

plot(corrWithMean)
%pause()

% Create an average over each frame, save as stackMean
stackSum2 = sum( stack(:,:,maxWindowCorrIndex:maxWindowCorrIndex+windowSize), 3);
stackMean2 = uint16(stackSum2/windowSize);
% Creates an array that holds correlation coefficients between stackMean
% and each frame of the stack
corrWithMean2 = zeros(frames,1);
for i = 1:frames
    corrWithMean2(i) = corr2(stack(:,:,i),stackMean2);
end

hold on
plot(corrWithMean2)
%pause()

[val ind] = sort(corrWithMean,'descend');
val(1:10);
ind(1:10);


calmTimeSt=maxWindowCorrIndex;
calmTimeEnd=maxWindowCorrIndex+windowSize;
M_calm = ind(1:10);
M_mov = ind( length(corrWithMean)-10: length(corrWithMean));

end

