
addpath('./');

datapoints = 100;
tLim = 10; % (In Seconds)
step = 1/datapoints;
fps = datapoints/tLim;

t = 0:step:tLim-step;
x = sin(2*pi*2*t) + sin(2*pi*4*t) + sin(2*pi*6*t);

low_pass = false;

Fourier(x, tLim, low_pass, datapoints)
