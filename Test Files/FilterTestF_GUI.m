function [ ] = FilterTestF_GUI( )
%FILTERTESTF_GUI Summary of this function goes here
%   Detailed explanation goes here
plot1 = false; % Data and Hanning window
plot2 = false; % Fourier transformed data
plot3 = false; % Sigmoid applied to FTdata
plot4 = true; % Filtered data
high_pass = false;
close all

slmin=1;slmax=30;
sliderVala=5;
sliderValb=5;
f = figure;
% pos(1) and pos(2) are the coordinate of the legend windows
% pos(3) and pos(4) are width and height of this windows
uicontrol('Style','slider','Callback',@sliderCallbacka,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderVala,...
                    'Position',[50 5 200 20]);
uicontrol('Style','slider','Callback',@sliderCallbackb,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderValb,...
                    'Position',[320 5 200 20]);
uicontrol('Button','button','Callback',@buttonCallback,'Position',[270 10 35 20]);
a = 10^(sliderVala/5);
b = 8*sliderValb;

sigmoidMir = 0;
                
function sliderCallbacka(src, evt)
    sliderVala=round(get(src, 'Value'));
    a = 10^(sliderVala/5)
    parameterTweaking();
end

function sliderCallbackb(src, evt)
    sliderValb=round(get(src, 'Value'));
    b = 10*sliderValb
    parameterTweaking();
end

function buttonCallback(src, evt)
    close all
    inverseTransform()
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
% time space "t" -> freqeuncy space "f'
f = (0:length(ft_x)-1)*(1/step)/length(ft_x);
f_leng = length(f);
% Plot fourier transformed data
if plot2
    plot(f,ft_x);
    hold on
    plot(f,real(ft_x));
    legend('x mirrored F.T.','real(x mirrored F.T.)')
    hold off
    pause()
end

fSig = f(1:length(f)/2);
parameterTweaking()

function parameterTweaking()
    uicontrol('Style','slider','Callback',@sliderCallbacka,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderVala,...
                    'Position',[50 5 200 20]);
    uicontrol('Style','slider','Callback',@sliderCallbackb,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderValb,...
                    'Position',[320 5 200 20]);
    uicontrol('Button','button','Callback',@buttonCallback,'Position',[270 10 35 20]);
    
    % The Sigmoid has the form (1+a e^(-f))^-b
    sigmoid = 1./(1+a*exp(-fSig)).^b;
    % Clever little way to find closest point to 0.5 in Sigmoid
    [delta index] = min(abs(0.5-sigmoid)); % Method 1
    % Closed form solution for cutoff Frequency
    %  High pass is flipped so cutoff will come from the right instead of
    %  from the left.
    f_cutoff = abs( log( (2^(1/b)-1)/a ) ); % Method 2

    if high_pass
        sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass
        cutoff_Freq1 = f(length(sigmoid)-index) % From indexing
        cutoff_Freq2 = max(fSig)-f_cutoff; % From closed form
    else
        sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
        cutoff_Freq1 = f(index);
        cutoff_Freq2 = f_cutoff;
    end
    cutoff_Freq = cutoff_Freq2;
    
    % Since we are dealing with mirrored data, we only plot the first half,
    %  the other half is redundant for GUI use
    half_f = f(1:f_leng/2);
    half_sigmoidMir = sigmoidMir(1:f_leng/2);
    half_ft_x = ft_x(1:f_leng/2);
    
    % Plot Sigmoid related things
    plot(half_f,real( half_sigmoidMir.*(max(half_ft_x)) )  );
    hold on
    plot(half_f,real( half_ft_x.*half_sigmoidMir )   );
    
    % Plot Vertical black line at cutoff frequency
    plot([cutoff_Freq cutoff_Freq], [-100 150],'linewidth',1.1,'color','k')
    % Add legend and plot labels
    legend('Normalized Sigmoid (*factor)','Sigmoid multiplied to F.T.x')
    title('Fourier Transformed Function with Hamming Window')
    xlabel('Frequency [Hz]') % x-axis label
    ylabel('Magnitude of Intensity') % y-axis label
    hold off
end

function inverseTransform()
    % Inverse Fourier Transform to get back filtered data
    yMir = ifft(ft_x.*sigmoidMir);
    y = yMir(tlen:2*tlen-1);
    if true
        plot(t,real(y));
        %pause(1)
        hold on
        plot(t,x,'linewidth',0.05);
        legend('filtered data','original data')
        pause()
    end
end

%close all

end