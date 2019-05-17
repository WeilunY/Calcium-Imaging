
addpath('./');

data = csvread('sample_data.csv');
s = size(data);
fps = 3.91;

datapoints = s(2);

low_pass = false;

Fourier2(data, fps, low_pass, datapoints)



