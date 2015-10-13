brain = im2double(imread('brain empty.png'));
brain_gray = im2double(rgb2gray(brain));
brain_bw = im2bw(brain);

cc = bwconncomp(brain_bw,4);

cortex = ones(size(brain_bw));
cortex(cc.PixelIdxList{1}) = 0;
cortex(cc.PixelIdxList{2}) = 0;

wm = zeros(size(brain_bw));
wm(cc.PixelIdxList{2}) = 1;

se = strel('disk',3); 
wm_eroded = imerode(wm,se);

% cables = im2double(rgb2gray(imread('wires.jpg')));
% cables = imresize(cables,size(brain_bw));
% cables(~wm) = 1;

cables = im2double(imread('wires.jpg'));
cables = imresize(cables,size(brain_bw));
cables = mask_rgb(cables,wm_eroded);

brain = color_mask(cortex,0.15,0.15,0.78);

background = imcomplement(wm + cortex - wm + wm_eroded);
background_color = color_mask(background,1,1,1);

imshow(background_color + brain + cables)

