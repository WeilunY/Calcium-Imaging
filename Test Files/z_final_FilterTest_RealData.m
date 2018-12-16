function [  ] = z_final_FilterTest_RealData( data )
%FILTERTESTF_GUI Summary of this function goes here
%   Detailed explanation goes here
plot1 = false; % Data and Hanning window
plot2 = false; % Fourier transformed data
plot3 = false; % Sigmoid applied to FTdata
plot4 = false;  % Filtered data
high_pass = true;
close all

% datapoints = 50;
% tLim = 10; % (In Seconds)
% step = 1/datapoints;
% t = 0:step:tLim-step;
% x = sin(2*pi*8*t) + sin(2*pi*22*t);

datapoints = size(data,2);
fps = 3.91;
tLim = datapoints/fps; % (In Seconds)
step = 1/fps;
t = 0:step:tLim-step;
x = data;

disp(['tLim: ',num2str(tLim)])
disp(['datapoints: ',num2str(datapoints)])
disp(['step: ',num2str(step)])

tlen = length(t);
disp(['tlen: ',num2str(tlen)])
mirLen = 3*tlen;
% Mirror your function to emulate periodicity
tMir = 0:step:(3*tLim)-step;
xMir = [ -fliplr(x), x, -fliplr(x) ];

slmin=0;slmax=50;
sliderVala=5;
sliderValb=5;
fig = figure('KeyPressFcn',@keypress,'units','pixels',...
              'position',[200 200 900 600],...
              'menubar','none',...
              'name','Fourier GUI',...
              'resize','off');
          % [DistFromBottom DistFromLeft PercentWidth PercentHeight]
ax1 = axes('Position',[0.1 0.2 0.8 0.7]);
% pos(1) and pos(2) are the coordinate of the legend windows
% pos(3) and pos(4) are width and height of this windows

a = 10;
b = 10^(sliderValb/20);

sigmoidMir = 0;
frequency_manually_set = false;
cutoff_f = 1;
                
% function sliderCallbacka(src, evt)
%     sliderVala=round(get(src, 'Value'));
%     a = 10^(sliderVala/5)
%     parameterTweaking();
% end

function sliderCallbackb(src, evt)
    sliderValb=round(get(src, 'Value'));
    b = 10^(sliderValb/20)
    parameterTweaking();
end

function buttonCallback(src, evt)
    close all
    inverseTransform()
end

function [] = textboxCallback(src, evt)
    text = get(src,'string');
    cutoff_f = str2double(text);
    
    if high_pass
        disp('CUTOFF')
        %disp( log2(a*exp((max(fSig)-cutoff_f))+1)^-1 )
        b = log2( a * exp( (cutoff_f-max(fSig)) )+1 )^-1;
    else
        disp('CUTOFF')
        disp( log2(a*exp(cutoff_f)+1)^-1 )
        b = log2(a*exp(cutoff_f)+1)^-1;
    end
    
    disp(['The cutoff frequency is: ',text])
    disp(['Parameter b set to: ',num2str(b)])
    frequency_manually_set = true;
    parameterTweaking();
end

% Apply a hanning window to your mirrored data
% This hanning window cuts off data a lot quicker
disp(['xMir size: ',num2str(size(xMir))])
disp(['mirLen: ',num2str(mirLen)])
disp(['size hann: ',num2str(size(transpose(hann(mirLen))))])

xMirHan = xMir.*transpose(hann(mirLen));
% This hanning window is a lot wider
big_hann = transpose(hann(3*mirLen));
size(big_hann(mirLen:2*mirLen))
size(xMir)
xMirHan = xMir.*big_hann(mirLen:2*mirLen-1);
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
% disp(ft_x)
% 
% plot(real(ft_x(1:50)))
% pause()


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
%     uicontrol(fig,'Style','slider','Callback',@sliderCallbacka,'Min',slmin,'Max',slmax,...
%                     'SliderStep',[1 1]./(slmax-slmin),'Value',sliderVala,...
%                     'Position',[50 5 200 20]);
    uicontrol(fig,'Style','slider','Callback',@sliderCallbackb,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderValb,...
                    'Position',[330 5 200 20]);
    uicontrol(fig,'Button','button','Callback',@buttonCallback,'Position',[270 10 45 30],...
                    'String','DONE','Tag','Button1');
    uicontrol(fig, 'Style','edit','Units','normalized','Position',[.8 .02 .06 .04],... 
                    'CallBack',@textboxCallback,'String','1');
    txt = uicontrol('Style','text','Units','normalized',...
                    'Position',[.7 .02 .1 .04],'String','Cutoff Frequency:');
    
    % The Sigmoid has the form (1+a e^(-f))^-b
    sigmoid = 1./(1+a*exp(-fSig)).^b;
    for i=1:size(fSig,2)
        ii = fSig(i);
        disp(ii)
        sigmoid(i) = heaviside(cutoff_f-ii);
        disp(sigmoid(i))
    end
    
    if frequency_manually_set
        cutoff_Freq = cutoff_f;
        
        if high_pass
            sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass
        else
            sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
        end
    else
        % Clever little way to find closest point to 0.5 in Sigmoid
        [delta index] = min(abs(0.5-sigmoid)); % Method 1
        % Closed form solution for cutoff Frequency
        %  High pass is flipped so cutoff will come from the right instead of
        %  from the left.
        f_cutoff = abs( log( (2^(1/b)-1)/a ) ); % Method 2

        if high_pass
            sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass
%             cutoff_Freq1 = f(length(sigmoid)-index) % From indexing
%             cutoff_Freq2 = max(fSig)-f_cutoff; % From closed form
            cutoff_Freq2 = abs(max(fSig)-f_cutoff); % From closed form
        else
            sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
%             cutoff_Freq1 = f(index);
            cutoff_Freq2 = f_cutoff;
        end
        cutoff_Freq = cutoff_Freq2;
    end
    
    % Since we are dealing with mirrored data, we only plot the first half,
    %  the other half is redundant for GUI use
    half_f = f(1:f_leng/2);
    half_sigmoidMir = sigmoidMir(1:f_leng/2);
    half_ft_x = ft_x(1:f_leng/2);
    
    disp(['f_leng: ',num2str(f_leng)])
    disp(['half_f size: ',num2str(size(half_f))])
    disp(['half_sigmoidMir size: ',num2str(size(half_sigmoidMir))])
    
    %plot( half_sigmoidMir.*abs(max(half_ft_x)) )
    %pause()
    
    % Plot Sigmoid related things
    plot(half_f,real( half_ft_x.*half_sigmoidMir )   );
    hold on
    plot(half_f,real( half_sigmoidMir.*abs(max(half_ft_x)) )  );
    
    
     
    % Plot Sigmoid related things X AXIS IS LOG
    %semilogx(half_f,real( half_sigmoidMir.*(max(half_ft_x)) )  );
    %hold on
    %semilogx(half_f,real( half_ft_x.*half_sigmoidMir )   );
    
    % Plot Vertical black line at cutoff frequency
    plot([cutoff_Freq cutoff_Freq], [-5 5],'linewidth',1.1,'color','k')
    % Add legend and plot labels
    legend('Normalized Sigmoid (*factor)','Sigmoid multiplied to F.T.x')
    title('Fourier Transformed Function with Hamming Window')
    xlabel('Frequency [Hz]') % x-axis label
    ylabel('Magnitude of Intensity') % y-axis label
    hold off
    frequency_manually_set = false;
end

function inverseTransform()
    ft_x = ft_x(1:size(sigmoidMir,2));
    
    disp(['ft_x size: ',num2str(size(ft_x))])
    disp(['sigmoidMir size: ',num2str(size(sigmoidMir))])
    
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