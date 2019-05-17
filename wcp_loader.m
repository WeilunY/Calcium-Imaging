function [ ] = wcp_loader( )
%TEST Summary of this function goes here
%   Detailed explanation goes here
root_wcp = '/Users/weilunyao/Desktop/Lab/WCP/';
root_mc = '/Users/weilunyao/Desktop/Lab/Pilo 10_6_17_cage3364992_mouse0/Data/';

%% Load WCP file contents
wcp_fp = strcat(root_wcp,'20171108_PiloMouse0_001[1-1].mat');
wcp_contents = load( wcp_fp );

wcp_t = wcp_contents.T1; % fps ~ 112.64 fps
wcp_t_len = length(wcp_t);
total_time_wcp = wcp_t(wcp_t_len);
wcp_y = wcp_contents.Y1;
wcp_y1 = wcp_contents.Y1(:,1);
wcp_y2 = wcp_contents.Y1(:,2);
wcp_y3 = wcp_contents.Y1(:,3);

wcp_t_spacing = wcp_t(2)-wcp_t(1);


disp('__WCP FILE__')
fprintf('WCP file is                : %d seconds long\n', total_time_wcp);
fprintf('WCP file fps is            : %d fps\n', wcp_t_len/total_time_wcp);
fprintf('WCP time file has length   : %d \n', wcp_t_len);
% fprintf('WCP time file has speed of : %f seconds per frame.\n', wcp_t_spacing);
fprintf('\n')

%% Load Motion-Cut output file
motion_cut_fp = strcat(root_mc,'20171108_PiloMouse0_001 motionCut_Output.mat');
motion_cut_contents = load(motion_cut_fp);

frameCorrStatus = motion_cut_contents.frameCorrStatus;% 3.91 fps
frameCorrStatus_frames = length(frameCorrStatus);
fps_movie = 3.91;
len_movie = frameCorrStatus_frames / fps_movie;

disp('__MOTION CUT FILE__')
fprintf('Cut movie has: %d frames.\n', frameCorrStatus_frames);
fprintf('Cut movie has: %d seconds.\n', len_movie);
fprintf('\n')

%% Cut Movie
Big_counter = 0;
small_counter = 0;
cut_data = [];
for i = frameCorrStatus
    Big_counter = Big_counter + 1;
    if Big_counter == 6
        while small_counter < 28
            cut_data = [cut_data, i];
            small_counter = small_counter + 1;
        end
    else
    
    while small_counter < 29
        cut_data = [cut_data, i];
        small_counter = small_counter + 1;
    end
    
    end
    small_counter = 0;
    
end

disp(length(cut_data))
        

cut_y1 = wcp_y1;
cut_y2 = wcp_y2;
cut_y3 = wcp_y3;

%% cut off data
for i = 1 : wcp_t_len
    if cut_data(i) == 0
        cut_y1(i) = 0;
        cut_y2(i) = 0;
        cut_y3(i) = 0;
    end
end

hold on 
plot( wcp_y1, 'r');
plot(cut_y1, 'b');

title('Y components in wcp');
xlabel('frame');
        
        
    

%% Analysis
wcp_frames_per_movie_frame = wcp_t_len/frameCorrStatus_frames;

disp('__ANALYSIS__')
fprintf('WCP frames per movie frame: %f \n', wcp_frames_per_movie_frame);
fprintf('')

%% Plotting
imagesc(frameCorrStatus)
%hold on
%plot(wcp_t,wcp_y1)
%hold off
end