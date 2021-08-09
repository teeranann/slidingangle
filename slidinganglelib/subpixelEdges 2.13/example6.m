%% SUBPIXEL EDGES - EXAMPLE 6 -----------
% SUBPIXEL EDGE DETECTION IN A REAL ANGIOGRAPHY

addpath(genpath('.'));

%% load image
url='http://serdis.dis.ulpgc.es/~atrujillo/ngImgCrop-master/test/angio2.PNG';
image = rgb2gray(imread(url));

%% subpixel detection
threshold = 4;
iter = 3;
[edges, RI] = subpixelEdges(image, threshold, 'SmoothingIter', iter); 

%% show image
showRestoredImage = true;
if showRestoredImage
    imshow(RI/255,'InitialMagnification', 'fit');
else
    imshow(image,'InitialMagnification', 'fit');
end

%% show edges
visEdges(edges);

