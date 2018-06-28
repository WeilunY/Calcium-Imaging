function [  ] = FilterTest1(  )
%FILTERTEST Summary of this function goes here
%   Detailed explanation goes here

W0 = 30;

%%% Generating the function
x = 0:pi/1000:20*pi;
l = length(x);
y = zeros(1,l);
for i=1:l
    xCor = i * pi/1000;
    w = (i^0.7)/10;
    y(i) = sin(w*x(i));
    if i==1 || i==l
        w
    end
end
conv = l/x(l);
disp(['initial angular frequency: ' , num2str( ((1)^0.7)/10 ) ]);
disp(['time=2 angular frequency: '  , num2str( ((2*conv)^0.7)/10 ) ]);
disp(['time=10 angular frequency: ' , num2str( ((10*conv)^0.7)/10 ) ]);
disp(['time=40 angular frequency: ' , num2str( ((40*conv)^0.7)/10 ) ]);
disp(['time=50 angular frequency: ' , num2str( ((50*conv)^0.7)/10 ) ]);
disp(['time=60 angular frequency: ' , num2str( ((60*conv)^0.7)/10 ) ]);
disp(['final angular frequency: '   , num2str( (l^0.7)/10 ) ]);

%%% Actual filterings
% Filters out signals above/below 1/tau
T = pi/1000; 
tau = 1/W0; %tau = 1/30;
a = T/tau;
% Where a = T/tau, T = the time between samples, and tau is the filter time constant.
yfiltHP = filter([1-a a-1],[1 a-1], y);
yfiltLP = filter(a, [1 a-1], y);

windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
% Find the moving average of the data and plot it against the original data.
y = filter(b,a,y);

%%% Plotting
hold off
plot(x,y,'b','linewidth',.5), grid on
hold on
plot(x,yfiltHP-2,'r','linewidth',.5)
plot(x,yfiltLP+2,'r','linewidth',.5)
%plot(x,y+5,'k','linewidth',.5)

end