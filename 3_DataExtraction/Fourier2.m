
% x: function to be transformed
% tLim: the time of the vieo
% low_pass: Boolean, whether user low or high pasw filter
% datapoints: the num of datapoints in video
function [] = Fourier(x, tLim, low_pass, datapoints)

plot1 = true; % Data and Hanning window
plot2 = true; % Fourier transformed data
plot3 = true; % Sigmoid applied to FTdata
plot4 = true;  % Filtered data

step = 1;
fps = datapoints/tLim;

t = 1:step:datapoints;
tlen = length(t);
mirLen = 3 * tlen;

%% Mirror function to emulate periodicity 
tMir = 0:step:(3*datapoints)-step;
xMir = [ fliplr(x), x, fliplr(x) ];

% GUI
slmin=0;slmax=50;
sliderValb=5;
fig = figure('KeyPressFcn',@keypress,'units','pixels',...
              'position',[200 200 900 600],...
              'menubar','none',...
              'name','Fourier GUI',...
              'resize','off');
          
        
a = 10;
b = 10^(sliderValb/20);

sigmoidMir = 0;
frequency_manually_set = false;
cutoff_f = 1;

%% GUI Callbacks

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
    
    if low_pass
        disp( log2(a*exp((max(fSig)-cutoff_f))+1)^-1 )
        b = log2( a * exp( (cutoff_f-max(fSig)) )+1 )^-1;
    else % High Pass
        disp( log2(a*exp(-1*cutoff_f)+1)^-1 )
        b = log2(a*exp(-1*cutoff_f)+1)^-1;
    end
    
    %disp(['The cutoff frequency is: ',text])
    %disp(['Parameter b set to: ',num2str(b)])
    frequency_manually_set = true;
    parameterTweaking();
end

%% Hanning Window
% Apply a hanning window to your mirrored data
% This hanning window cuts off data a lot quicker
%xMirHan = xMir.*transpose(hann(mirLen));


% This hanning window is a lot wider
big_hann = transpose(hann(3*mirLen));
size(big_hann(mirLen:2*mirLen))
size(xMir)
xMirHan = xMir.*big_hann(mirLen:2*mirLen-1);
if plot1
    plot(tMir,xMir)
    hold on
    plot(tMir,xMirHan)
    plot([datapoints datapoints], [-1.5 1.5],'linewidth',1.5,'color','k')
    plot([2*datapoints 2*datapoints], [-1.5 1.5],'linewidth',1.5,'color','k')
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


% Power specturm
ft_x_2 = ft_x(1:f_leng/2) .* conj(ft_x(1:f_leng/2));
ft_x_s = sum(ft_x_2);
power_spectrum = ft_x_2 / ft_x_s;

if true
    plot(f(1:f_leng/2), power_spectrum);
    hold on
    plot(f(1:f_leng/2), power_spectrum);
    title('Normalized Power Spectrum');
    xlabel('Frequency');
    ylabel('Percentage');
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
    
    
    if frequency_manually_set
        cutoff_Freq = cutoff_f;
        
        if low_pass
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

        if low_pass
            sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass
            %cutoff_Freq1 = f(length(sigmoid)-index) % From indexing
            cutoff_Freq2 = max(fSig)-f_cutoff; % From closed form
            cutoff_Freq2 = abs(max(fSig)-f_cutoff); % From closed form
        else
            sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
            %cutoff_Freq1 = f(index); % From indexing
            cutoff_Freq2 = f_cutoff; % From closed form
            disp('CUTOFF FREQ (PARAMETER TWEAKING)')
            disp(cutoff_Freq2)
        end
        cutoff_Freq = cutoff_Freq2;
    end
    
    % Since we are dealing with mirrored data, we only plot the first half,
    %  the other half is redundant for GUI use
    
    %% Note: un do half
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
    frequency_manually_set = false;
end

function inverseTransform()
    % Inverse Fourier Transform to get back filtered data
    yMir = ifft(ft_x.*sigmoidMir);
    yMir = yMir./ big_hann(mirLen:2*mirLen-1);
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