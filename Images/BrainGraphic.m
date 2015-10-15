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
% cortex_filtered = imgaussfilt(cortex_mask,2);
cortex_filtered = imfilter(cortex_mask,fspecial('gaussian'));
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

cables_edge = color_mask(cables_edge_mask,170/255,56/255,50/255);


brain = color_mask(cortex_mask,198/255,56/255,32/255);

% subplot(1,2,1)
% imshow(brain + background_color + repmat(wm,[1 1 3]))
% subplot(1,2,2)
% imshow(brain + cables_edge + cables_filling + background_color)

imwrite(brain + cables_edge + cables_filling + background_color,'cables.png')
imwrite(brain + background_color + repmat(wm,[1 1 3]),'cortex.png')

N=40;
n=3;
nn=1;

box = (-1)*ones(N+2*n);
box(nn+1:N+2*n-nn,nn+1:N+2*n-nn) = 1;
box(n+1:N+n,n+1:N+n) = 0;

mask = cortex_bw + wm;
grid = zeros(size(cortex_bw));
i=n+1;
while i+N-1 < size(cortex_bw,1)
   j=n+1;
   while j+N-1 < size(cortex_bw,2)
    if sum(sum(mask(i:i+N-1,j:j+N-1))) > 0
       grid(i-n:i+N+n-1,j-n:j+N+n-1) = box;
       
       if sum(sum(mask(i-N-2*n:i-2*n-1,j:j+N-1))) == 0
           grid(i-n-1,j-n:j+N+n-1) = -1;
       end
       
       if sum(sum(mask(i+N+2*n:i+2*N+2*n-1,j:j+N-1))) == 0
           grid(i+N+1,j-n:j+N+n-1) = -1;
       end
       
       if sum(sum(mask(i:i+N-1,j+N+2*n:j+2*N+2*n-1))) == 0
       grid(i-n-1:i+N+n-1,j+N+n) = -1;
       end
       
       if sum(sum(mask(i:i+N-1,j-N-2*n:j-2*n-1))) == 0
       grid(i-n-1:i+N+n-1,j-n-1) = -1;
       end
       
    end
    j = j + N + 2*n;   
   end
    i = i + N + 2*n;
end

% t=10;
% grid = [grid(:,t:end) grid(:,1:t-1)];

grid = repmat(grid,[1 1 3]);

imwrite(brain + cables_edge + cables_filling + background_color + grid,'grid.png')

%%
im=cables_edge_mask(200:250,200:250);

subplot(1,2,1)
imshow(im)
subplot(1,2,2)
imshow(im2ell(im))

%%

img = im2double(cortex_bw + cables_edge_mask);

ell_grid = zeros(size(cortex_bw));

ellipses = cell(1,9);

t = linspace(0,2*pi);
ell_coo = [cos(t)', sin(t)'];
ellipses{1} = coo2im(ell_coo,N,N);
imshow(ellipses{1})

a=1.5;b=1;
for d = 0:7
    phi = d*pi/8;
    ell_coo = [a*cos(phi)*cos(t)' - b*sin(phi)*sin(t)', a*sin(phi)*cos(t)' + b*cos(phi)*sin(t)'];
    subplot(3,3,d+1)
%     plot(ell_coo(:,1),ell_coo(:,2))
%     axis([-3 3 -3 3])
    ellipses{d+2} = coo2im(ell_coo,N,N);
    imshow(ellipses{d+2})
end

% 
% i=n+1;
% while i+N-1 < size(cortex_bw,1)
%    j=n+1;
%    while j+N-1 < size(cortex_bw,2)
%     if sum(sum(double(cortex_bw(i:i+N-1,j:j+N-1) + cables_edge_mask(i:i+N-1,j:j+N-1)))) > 0
%         
%         dir = zeros(1,4);
%         for d = 0:3
%         img_rot = imrotate(img,d*pi/4);
%         block = img_rot(i:i+N-1,j:j+N-1);
%         dir(1,i) = sum(sum(edge(block,'Prewitt',0.1,'horizontal')));
%         end
%         dir = dir/sum(dir);
% %        ell_grid(i:i+N-1,j:j+N-1) = im2ell(im2double(cortex_bw(i:i+N-1,j:j+N-1) + cables_edge_mask(i:i+N-1,j:j+N-1)));   
%      
%        
%     end
%     j = j + N + 2*n;   
%    end
%     i = i + N + 2*n;
% end
% 
% % ell_grid = [ell_grid(:,t:end) ell_grid(:,1:t-1)];
% 
% imshow(repmat(-0.5*ell_grid,[1 1 3]) + brain + cables_edge + cables_filling + background_color + grid)













