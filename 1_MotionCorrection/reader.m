%%%%just few lines of code to read your multitiff type of data
function [stack]=reader(fname)
info = imfinfo(fname);
num_images = numel(info);

imgTMP = round(imread(fname, 2)); % Temporarily holds the stack

stack = zeros(size(imgTMP,1), size(imgTMP,2), round(num_images/2), 'uint16');
h = waitbar(0,'Loading images.');% Generates loading bar

for index=2:2:num_images % Every other frame
%for index=1:1:num_images %%% If not every other frame
    waitbar(index/ num_images)
    imgTMP = imread(fname, index);
    imgTMP = uint16(imgTMP*64); % From 10bit to 16bit
    stack(:,:,index/2)=imgTMP; % Takes every OTHER frame
%    stack(:,:,index)=imgTMP; %%% If not every other frame
end
close(h)


end