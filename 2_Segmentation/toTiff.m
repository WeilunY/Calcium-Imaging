function[] = toTiff(stack,name)
%Converts a 3D matrix into a tiff stack

outputFileName = [name,'.tif'];
img = stack; 
frames = length(img(1, 1, :));

h = waitbar(0,'Converting Data into a CUT tiff file');

for K=1:frames
   imwrite(img(:, :, K), outputFileName, 'WriteMode', 'append');
   waitbar(K/frames)
end
close(h)

end

% toTiff(stackAdjustedCut,'name')