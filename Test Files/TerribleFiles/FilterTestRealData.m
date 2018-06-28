function [ ] = FilterTestRealData( Data )
%FILTERTEST2 FilterTestRealData(intensityData(n,:))

dataPoints = length(Data);
t = linspace(0,dataPoints-1,dataPoints);
x = Data;
%plot(x)
%pause()

tMir = linspace(-dataPoints,2*dataPoints,300);
xMir = [-fliplr(x), x, -fliplr(x)];
%xMir = [x, x, x];

plot(xMir)
pause()

% Real sine has f of about 0.038
lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',0.2,'PassbandRipple',.001, ...
         'SampleRate',3.81);
%fvtool(lpFilt)
yLP = filter(lpFilt,xMir);

% Noise has f of about 0.38
hpFilt = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',0.0065,'PassbandRipple',.04, ...
         'SampleRate',3.81);
%fvtool(hpFilt)
yHP = filter(hpFilt,xMir);

% Keeps frequencies between cutoff1 and cutoff2
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',0.001,'CutoffFrequency2',0.05, ...
         'SampleRate',3.81);
%fvtool(bpFilt)
yBP = filter(bpFilt,xMir);

length(yLP)*1/3;
length(yLP)*2/3;
yLP = yLP(length(yLP)*1/3:length(yLP)*2/3-1);
yHP = yHP(length(yHP)*1/3:length(yHP)*2/3-1);
yBP = yBP(length(yBP)*1/3:length(yBP)*2/3-1);

max(yHP)
max(x)

% Finds the maximizing values of the filters
%{l
h = waitbar(0,'Please wait...');
steps = 1000;
maxes = ones(steps,1);
for i=1:steps
    hpFilt = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',0.1/i,'PassbandRipple',0.006, ...
         'SampleRate',3.81);
    yHP = filter(hpFilt,xMir);
    maxes(i) = max(yHP);
    waitbar(i / steps);
end
close(h) 
plot(maxes)
max(maxes)
pause();
%}

%Plot data
hold off
plot(t,x)
hold on
%plot(t,yLP+0.5);
plot(t,yHP+1);
plot(t,x-yHP+2);
%plot(t,yBP+1.5)
%legend('Input','Low Pass','High Pass','Band Pass')
legend('Input','High Pass','difference')

end

