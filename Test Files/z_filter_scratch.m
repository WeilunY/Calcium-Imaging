function [  ] = z_filter_scratch( data )

fps = 3.91;
datapoints = size(data,2);
fft_data = fft(data);
filtered_fft_data = fft_data; % TEMPORARY

Fs = fps;
L = datapoints;
f = Fs*(0:(L))/L;


disp(fft_data)

freq_cutoff = 1;
slmin=0;slmax=50;
sliderVal = 10;
fig = figure('KeyPressFcn',@keypress,'units','pixels',...
              'position',[200 200 900 600],...
              'menubar','none',...
              'name','Fourier GUI',...
              'resize','off');
uicontrol(fig,'Style','slider','Callback',@sliderCallback,'Min',slmin,'Max',slmax,...
                    'SliderStep',[1 1]./(slmax-slmin),'Value',sliderVal,...
                    'Position',[330 5 200 20]);
uicontrol(fig,'Button','button','Callback',@buttonCallback,'Position',[270 10 45 30],...
                    'String','DONE','Tag','Button1');
                
plot_data()
                
function sliderCallback(src, evt)
    sliderVal=round(get(src, 'Value'));
    freq_cutoff = sliderVal/10;
    plot_data()
end

function buttonCallback(src, evt)
    final_plot()
%     close all
end


function plot_data()
    
%     plot(data)
%     pause(1)
    
    step_fcn = f;
    for i=1:size(f,2)
        ii = f(i);
        step_fcn(i) = heaviside(ii-freq_cutoff);
    end
    
    disp(size(step_fcn))
    disp(size(fft_data))
    step_fcn = step_fcn(1:size(fft_data,2));
    
    filtered_fft_data = fft_data.*step_fcn;
    
    plot( real(filtered_fft_data) )
%     pause()
    
end

function final_plot()
    
    filtered_data = ifft(filtered_fft_data);
    
    plot( data )
    hold on
    plot( real(filtered_data) )
    hold off
    legend('DATA','FILTERED DATA')
    
end



end

