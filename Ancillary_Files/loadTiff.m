function [  ] = loadTiff( filePath, fileNameTiff, fileName  )
%LOADTIFF Summary of this function goes here
%   Detailed explanation goes here
filePathSaveTo = (strcat(filePath, 'Data/', fileName, ' stack'))
filePathToTiff = strcat(filePath,fileNameTiff);

tiff_info = imfinfo(filePathToTiff); % return tiff structure, one element per image
tiff_stack = imread(filePathToTiff, 1) ; % read in first image
%concatenate each successive tiff to tiff_stack
for ii = 2 : size(tiff_info, 1)
    temp_tiff = imread(filePathToTiff, ii);
    tiff_stack = cat(3 , tiff_stack, temp_tiff);
end
stack = tiff_stack;

save(char(filePathSaveTo),'stack'); % Saves into .mat file ending in 'stack'

end
