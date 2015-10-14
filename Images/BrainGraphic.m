brain = im2double(imread('brain empty.png'));
brain_gray = im2double(rgb2gray(brain));
brain_bw = im2bw(brain);

cc = bwconncomp(brain_bw,4);

cortex_mask = ones(size(brain_bw));
cortex_mask(cc.PixelIdxList{1}) = 0;
cortex_mask(cc.PixelIdxList{2}) = 0;
cortex_mask = imdilate(cortex_mask,strel('disk',2));

cortex_edge = edge(cortex_mask,'Canny');
cortex_edge_dilated = imdilate(cortex_edge,strel('disk',2));
cortex_filtered = imgaussfilt(cortex_mask,2);
cortex_mask(cortex_edge_dilated) = cortex_filtered(cortex_edge_dilated);

cortex_bw = cortex_mask > 0;

cc2 = bwconncomp(~cortex_bw,4);

wm = ones(size(brain_bw));
wm(cc2.PixelIdxList{1}) = 0;
wm(cortex_bw) = 0;
wm_eroded = imerode(wm,strel('disk',5));

background = ones(size(brain_bw));
background(cc2.PixelIdxList{2}) = 0;
background(cortex_bw) = 0;
background_color = color_mask(background,250/255,230/255,230/255);

cables = im2double(imread('Spaghetti.jpg'));
cables = imresize(cables,size(brain_bw));
cables = mask_rgb(cables,wm);

cables_edge_mask = imdilate(edge(im2double(rgb2gray(cables)),'canny'),strel('disk',1));
cables_edge_mask(logical(imcomplement(wm))) = 0;
cables_edge_mask(logical(imdilate(cortex_mask,strel('disk',3)))) = 0;

cables_filling = imcomplement(cables_edge_mask);
cables_filling(logical(background)) = 0;
cables_filling(logical(cortex_bw)) = 0;
cables_filling = color_mask(cables_filling,1,1,1);

cables_edge = color_mask(cables_edge_mask,0.25,0.25,0.25);


brain = color_mask(cortex_mask,198/255,56/255,32/255);

subplot(1,2,1)
imshow(brain + background_color + repmat(wm,[1 1 3]))
subplot(1,2,2)
imshow(brain + cables_edge + cables_filling + background_color)

imwrite(brain + cables_edge + cables_filling + background_color,'cables.png')
imwrite(brain + background_color + repmat(wm,[1 1 3]),'cortex.png')
