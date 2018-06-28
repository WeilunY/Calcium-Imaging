function [ ] = FilterTest3( )
%FILTERTEST2 Summary of this function goes here
%   Detailed explanation goes here

% Create a 1-by-100 row vector of sinusoidal data that is corrupted by random noise.
t = linspace(-pi,pi,100);
rng default  %initialize random number generator
noise = 0.25*rand(size(t));
%x = sin(t) + 0.25*rand(size(t));
x = sin(t) + noise;

tMir = linspace(-2*pi,2*pi,300);
xMir = [fliplr(x), x, fliplr(x)];
xMir = [x, x, x];

% Real sine has f of about 0.038
lpFilt = designfilt('lowpassiir','FilterOrder',8, ...
         'PassbandFrequency',0.2,'PassbandRipple',.001, ...
         'SampleRate',3.81);
fvtool(lpFilt)
yLP = filter(lpFilt,xMir);

% Noise has f of about 0.38
hpFilt = designfilt('highpassiir','FilterOrder',8, ...
         'PassbandFrequency',0.1,'PassbandRipple',0.5, ...
         'SampleRate',3.81);
fvtool(hpFilt)
yHP = filter(hpFilt,xMir);

% Keeps frequencies between cutoff1 and cutoff2
bpFilt = designfilt('bandpassfir','FilterOrder',20, ...
         'CutoffFrequency1',0.001,'CutoffFrequency2',0.05, ...
         'SampleRate',3.81);
fvtool(bpFilt)
yBP = filter(bpFilt,xMir);

length(yLP)*1/3;
length(yLP)*2/3;
yLP = yLP(length(yLP)*1/3:length(yLP)*2/3-1);
yHP = yHP(length(yHP)*1/3:length(yHP)*2/3-1);
yBP = yBP(length(yBP)*1/3:length(yBP)*2/3-1);

size(yHP)


x_fft = abs(fft(xMir)); 
plot(x_fft)
pause();


%Plot data
hold off
plot(t,x);
hold on
plot(t,yLP);
plot(t,yHP);
plot(t,yBP)
plot(t,noise);
legend('Input','Low Pass','High Pass','Band Pass','noise')

end

