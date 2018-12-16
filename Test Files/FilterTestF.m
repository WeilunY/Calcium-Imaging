function [ ] = FilterTestF( )
%FILTERTESTF Summary of this function goes here
%   Detailed explanation goes here
plot1 = true; % Data and Hanning window
plot2 = true; % Fourier transformed data
plot3 = true; % Sigmoid applied to FTdata
plot4 = true; % Filtered data
high_pass = true;

% Coefficents for Sigmoid function
if high_pass
    a = 500; b = 5; % HP
else
    a = 10000000; b = 10; % LP
end
step = 1/50;
tLim = 10;

t = 0:step:tLim-step;
x = sin(2*pi*15*t) + sin(2*pi*20*t);

tlen = length(t)
mirLen = 3*tlen;

% Mirror your function to emulate periodicity
tMir = 0:step:(3*tLim)-step;
xMir = [ -fliplr(x), x, -fliplr(x) ];

% Apply a hanning window to your mirrored data
xMirHan = xMir.*transpose(hann(mirLen));
if plot1
    plot(tMir,xMir)
    hold on
    plot(tMir,xMirHan)
    plot([tLim tLim], [-1.5 1.5],'linewidth',1.5,'color','k')
    plot([2*tLim 2*tLim], [-1.5 1.5],'linewidth',1.5,'color','k')
    legend('x Mirrored','x Mirrored Hanning')
    hold off
    pause()
end

% Apply fast fourier transform to the data
ft_x = fft(xMirHan);
f = (0:length(ft_x)-1)*(1/step)/length(ft_x); % t -> f
% Plot fourier transformed data
if plot2
    plot(f,ft_x);
    hold on
    plot(f,real(ft_x));
    legend('x mirrored F.T.','real(x mirrored F.T.)')
    hold off
    pause()
end

% Create sigmoid function
fSig = f(1:length(f)/2);
sigmoid = 1./(1+a*exp(-fSig)).^b;
if high_pass
    sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass
else
    sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
end
% Plot Sigmoid related things
if plot3
    plot(f,real( sigmoidMir.*(max(ft_x)) ));
    hold on
    plot(f,real( ft_x.*sigmoidMir ));
    legend('Normalized Sigmoid (*factor)','Sigmoid multiplied to F.T.x')
    hold off
    pause()
end

% Inverse Fourier Transform to get back filtered data
yMir = ifft(ft_x.*sigmoidMir);
y = yMir(tlen:2*tlen-1);
if plot4
    plot(t,real(y));
    %pause(1)
    hold on
    plot(t,x,'linewidth',0.05);
    legend('filtered data','original data')
    pause()
end

close all

end