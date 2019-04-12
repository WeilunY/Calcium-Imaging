
addpath('./');

data = csvread('sample_data.csv');
s = size(data);
fps = 3.91;

datapoints = s(2);
tLim = datapoints/fps; % (In Seconds)

step = 1/datapoints;


t = 0:step:tLim-step;
low_pass = false;

Fourier2(data, tLim, low_pass, datapoints)

