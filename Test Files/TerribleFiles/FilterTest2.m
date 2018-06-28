function [ ] = FilterTest2( )
%FILTERTEST2 Summary of this function goes here
%   Detailed explanation goes here

% Create a 1-by-100 row vector of sinusoidal data that is corrupted by random noise.
t = linspace(-pi,pi,100);
rng default  %initialize random number generator
noise = 0.25*rand(size(t));
%x = sin(t) + 0.25*rand(size(t));
x = sin(t) + noise;
% For a window size of 5, compute the numerator and denominator coefficients 
% for the rational transfer function.
windowSize = 5; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
% Find the moving average of the data and plot it against the original data.
yfiltLP = filter(b,a,x); % Should become the sin wave
yfiltHP = filter([0.8000 -0.8000],a,x); % Should become the noise
% err=1.5

%{s
lowestErr = 1.6;
goodCoef = [0 0];
ttt = zeros(1,1000);
for a1=1:1000
    for a2=-1000:0
        yfiltHP = filter([a1/1000 a2/1000],a,x);
        % Calculate error
        LPerror = 0;
        HPerror = 0;
        % Least squares error analysis
        HPerror = sum((yfiltHP-noise).^2)^0.5;
        % Records the best values we find
        if HPerror<lowestErr
            lowestErr=HPerror;
            goodCoef = [a1 a2];
        end
        if a1==425 && a2==-393
            ttt(a2+801) = HPerror;
        end
    end
end
lowestErr
goodCoef
%}
yfiltHP = filter(goodCoef/1000,a,x); % Should become the noise


%Plot data
hold off
plot(t,x)
%hold on
plot(t,yfiltLP)

plot(t,yfiltHP-mean(yfiltHP),'linewidth',.5)
hold on
plot(t,noise-mean(noise),'linewidth',.5)
%plot(t,yfiltHP-noise,'linewidth',.5)
legend('HP Filtered Data','pure noise')
%legend('Input Data','LP Filtered Data','HP Filtered Data','pure noise')

end

