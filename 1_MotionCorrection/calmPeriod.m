function [calmTimeSt,calmTimeEnd,expStarted,expEnded, Y1, M_calm, M_mov]=calmPeriod(movementSF, dataFile, imdilateParameter)
load(dataFile);
%expStarted=T1(1)    %find(Y1(:,3)>mean(Y1(:,3)),1,'first')/movementSF;
%expEnded=T1(length(T1))     %find(Y1(:,3)>mean(Y1(:,3)),1,'last')/movementSF;
expStarted=find(Y1(:,3)>mean(Y1(:,3)),1,'first')/movementSF;
expEnded=find(Y1(:,3)>mean(Y1(:,3)),1,'last')/movementSF;


%%%%% Sampling frequency of the motion
seconds = round(expEnded-1);
secondPart=zeros(1,seconds); % Partitioning to seconds
for index=1:(seconds-1) % The last section can be incomplete
    secondPart(index)=std(Y1(movementSF*(index-1)+1:movementSF*index,1));
end

% I = find(secondPart==0,1,'first'); %%%%For some unknown reason there might be some zeros at the end of the data
% secondPart(I:end)=[];

stdVect=secondPart'; 
stdVect(:,2)=linspace(1,length(stdVect),length(stdVect))'; % 1-300 in steps of 1
motionSorted=sortrows(stdVect);%%%% This line and the following one are faster than "prctile" 
                                        %%% function from Statistics Toolbox from Mathworks


M_calm=motionSorted(1:10,2);
M_mov=motionSorted(end-9:end,2);
motionSorted(:,2)=[];


THR=motionSorted(round(length(motionSorted)*0.2)); %%%We seek one tenth of the time; the calmest moments
motionOver=find(secondPart>THR);
secondPartBinary=secondPart;
secondPartBinary(motionOver)=0;

calmMomentsBinary=zeros(size(secondPartBinary));
calmMomentsBinary(find(secondPartBinary))=1;

imdilatePattern=ones(1,imdilateParameter); %%% Not everything is perfect on this world. 
%%% Here we fill points like "mouse makes  a tiny move (but over THR) during otherwise long calm period"
secondPartBinary5=imdilate(calmMomentsBinary,imdilatePattern); %%%Filling the gaps

CC = bwconncomp(secondPartBinary5);
[~, segmentN]=size(CC.PixelIdxList);
segmentsSizes=zeros(1,segmentN);
for index=1:segmentN %%%We measure the size of the identified calm periods
 segmentsSizes(index)=length(CC.PixelIdxList{1,index});   
end
[~,I] = max(segmentsSizes); %%%%%Which calm period is biggest
calmSeconds=CC.PixelIdxList{1,I};
calmTimeSt=(calmSeconds(1));
calmTimeEnd=(calmSeconds(end));

plot(T1, Y1(:,1))
grid minor
maxValue=max(Y1(:,1));
minValue=min(Y1(:,1));
x1=[calmTimeSt calmTimeSt+0.001];
x2=[calmTimeEnd calmTimeEnd+0.001];
hold on
plot(x1, [minValue, maxValue], 'r', 'LineWidth', 1.5);
plot(x2, [minValue, maxValue], 'r', 'LineWidth', 1.5);

%pause;
close;

end