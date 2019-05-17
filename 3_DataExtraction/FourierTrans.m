function [] = FourierTrans(x, fps, low_pass, datapoints)

%% Data setup:
step = 1/fps;
tLim = datapoints / fps;

t = 0: step : tLim - step;
tlen = length(t);
mirLen = 3 * tlen;

%% Mirrored Data:
tMir = 0:step:(3*tLim)-step;
xMir = [ fliplr(x), x, fliplr(x) ];


sigmoidMir = 0;
cutoff_f = 0.3;

a = 10;
b = 10 ^(0.3/20);

%% Hanning Window
big_hann = transpose(hann(3*mirLen));
size(big_hann(mirLen:2*mirLen))
size(xMir)
xMirHan = xMir.*big_hann(mirLen:2*mirLen-1);

% fast fourier
ft_x = fft(xMirHan);

% time space "t" -> freqeuncy space "f'
f = (0:length(ft_x)-1)*(1/step)/length(ft_x);
f_leng = length(f);


%% Power specturm
ft_x_2 = ft_x(1:f_leng/2) .* conj(ft_x(1:f_leng/2));
ft_x_s = sum(ft_x_2);
power_spectrum = ft_x_2 / ft_x_s;


fSig = f(1:length(f)/2);

%% tweaking
sigmoid = 1./(1+a*exp(-fSig)).^b;
f_cutoff = abs( log( (2^(1/b)-1)/a ) );

if low_pass
    sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass
else
    sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
end

if low_pass
        sigmoidMir = [ fliplr(sigmoid), sigmoid ]; % High Pass

        cutoff_Freq2 = abs(max(fSig)-f_cutoff); % From closed form
        disp('CUTOFF FREQ (PARAMETER TWEAKING)')
        disp(cutoff_Freq2)
    else
        sigmoidMir = [ sigmoid, fliplr(sigmoid) ]; % Low Pass
        cutoff_Freq2 = cutoff_f; % From closed form
        disp('CUTOFF FREQ (PARAMETER TWEAKING)')
        disp(cutoff_Freq2)
end

cutoff_f = cutoff_Freq2;
 
 %% Note: un do half
    half_f = f(1:f_leng/2);
    half_sigmoidMir = sigmoidMir(1:f_leng/2);
    half_ft_x = ft_x(1:f_leng/2);
    
    % Plot Sigmoid related things
    plot(half_f,real( half_sigmoidMir.*(max(half_ft_x)) )  );
    hold on
    plot(half_f,real( half_ft_x.*half_sigmoidMir )   );
    
    % Plot Vertical black line at cutoff frequency
    plot([cutoff_f cutoff_f], [-5, 5],'linewidth',1.1,'color','k')
    % Add legend and plot labels
    legend('Normalized Sigmoid (*factor)','Sigmoid multiplied to F.T.x')
    title('Fourier Transformed Function with Hanning Window')
    xlabel('Frequency [Hz]') % x-axis label
    ylabel('Magnitude of Intensity') % y-axis label
    hold off
    
 %% Inverse Transfotm
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