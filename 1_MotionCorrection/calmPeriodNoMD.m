function [ calmTimeSt,calmTimeEnd, M_calm, M_mov, reference ] = calmPeriodNoMD( stack, fps )
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
windowSize = 50; % The number of frames the reference is taken over
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
legend({'correlation with stack-avg','correlation with stable-stack-subset-avg'},'Location','southwest')
xlabel('Frame') 
ylabel('Correlation Coefficient') 
%pause()
close()

% sorts the correlation with the mean in descending order. val contains the
% values in descending order. ind lists the indices of these values in the
% same order.
[val ind] = sort(corrWithMean,'descend');
val(1:10);
ind(1:10);


calmTimeSt=maxWindowCorrIndex;
calmTimeEnd=maxWindowCorrIndex+windowSize;
%pause()
M_calm = ind(1:10);
M_mov = ind( length(corrWithMean)-10: length(corrWithMean));
M_calm(1)
M_calm = M_calm./fps;
M_mov = M_mov./fps;
M_calm(1)

reference_cumunulative = double(zeros( 128, 512));

% Calculating the reference image, take average over the calm period we
% decided on. From frame = maxWindowCorrIndex to frame = maxWindowCorrIndex + windowSiz
for frame = maxWindowCorrIndex : maxWindowCorrIndex+windowSize % Start with frame 2.
  % Add each frame together into reference_cumunulative
  reference_cumunulative = reference_cumunulative + double(stack(:,:,frame)); 
end
% At the end divide by the number of frames we added together (windowSize)
reference = reference_cumunulative./windowSize;


end

