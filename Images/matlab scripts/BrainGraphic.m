cables = im2double(imread('images/Spaghetti.jpg'));

brain = im2double(imread('images/brain empty.png'));
brain_gray = im2double(rgb2gray(brain));
brain_bw = im2bw(brain);

brain_gray = imresize(brain_gray,[size(cables,1), size(cables,2)] );
brain_bw = imresize(brain_bw,[size(cables,1), size(cables,2)]);

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
background_color = color_mask(background,241/255,220/255,219/255); % BACKGROUND COLOR 250/255,230/255,230/255

% cables = imresize(cables,size(brain_bw));
cables = mask_rgb(cables,wm);

cables = imresize(cables,[round(size(cables,1)/2) round(size(cables,2)/2)]);
cables = im2double(rgb2gray(cables));

cables_filt = imfilter(cables,fspecial('disk',6));
cables(cables>0.1) = cables_filt(cables>0.1);

cables_edge_mask = edge(cables,'canny',0.03);
cables_edge_mask = imdilate(cables_edge_mask,strel('disk',1));
cables_edge_mask = imresize(cables_edge_mask,size(brain_bw));

cables_edge_mask(logical(imcomplement(wm))) = 0;
cables_edge_mask(logical(imdilate(cortex_mask,strel('disk',10)))) = 0;

cables_filling = imcomplement(cables_edge_mask);
cables_filling(logical(background)) = 0;
cables_filling(logical(cortex_bw)) = 0;
cables_filling = color_mask(cables_filling,1,1,1);

cables_edge = color_mask(cables_edge_mask,218/255,150/255,149/255); % CABLE COLOR 170/255,56/255,50/255


brain = color_mask(cortex_mask,149/255,55/255,53/255); % BRAIN COLOR 198/255,56/255,32/255

% subplot(1,2,1)
% imshow(brain + background_color + repmat(wm,[1 1 3]))
% subplot(1,2,2)
% imshow(brain + cables_edge + cables_filling + background_color)

imwrite(brain + background_color + repmat(wm,[1 1 3]),'output/1_cortex.png')
imwrite(brain + cables_edge + cables_filling + background_color,'output/2_cables.png')

%% GRID
N=160;
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

imwrite(brain + cables_edge + cables_filling + background_color + grid,'output/3_grid.png')


%% ELLIPSE TEMPLATES

base_ellipse = cell(1,4);
for i=1:3
    base_ellipse{i} = zeros(N);
end

for e = 1:3  
    for i=1:N
        for j=1:N
            if abs((N-i+1 - N/2)^2/(40-(e-1)*12)^2 + (j-N/2)^2/(40+(e-1)*12)^2 - 1) <= 0.05 + (e-1)^1.5*0.05
                base_ellipse{e}(i,j)=1;
            end      
        end
    end
    
    base_ellipse{e} = imdilate(base_ellipse{e},strel('disk',1));
end

rot_ell = cell(1,2*8+1);
rot_ell{1} = base_ellipse{1};
for k=1:2
    for r=0:7   
        rot_ell{(k-1)*8 + r + 2} = imrotate(base_ellipse{k+1},r*180/8,'nearest','crop');
    end
end

ell_hist = cell(1,2*8 + 1);
for r = 1:2*8+1
    [~,gdir] = imgradient(rot_ell{r});
    ell_hist{r} = hist( reshape(gdir(gdir~=0),1,numel(gdir(gdir~=0))) , linspace(-180,180,10) );
    ell_hist{r} = ell_hist{r} / sum(ell_hist{r});
end

%%

% for e = 1:3
% base_ellipse{e} = imdilate(base_ellipse{e},strel('disk',2)) - 0.5*base_ellipse{e};
% end
% 
% rot_ell = cell(1,2*8+1);
% rot_ell{1} = base_ellipse{1};
% for k=1:2
%     for r=0:7   
%         rot_ell{(k-1)*8 + r + 2} = imrotate(base_ellipse{k+1},r*180/8,'nearest','crop');
%     end
% end
% 
% imwrite(base_ellipse{1},'base1.png')
% imwrite(base_ellipse{2},'base2.png')
% imwrite(base_ellipse{3},'base3.png')

%%

base_ellipse{1} = im2double(imread('tensors/base1.png'));
base_ellipse{2} = im2double(imread('tensors/base2.png'));
base_ellipse{3} = im2double(imread('tensors/base3.png'));
base_ellipse{4} = im2double(imread('tensors/base0.png'));

rot_ell = cell(1,2*8+1);
rot_ell{1} = base_ellipse{1};
for k=1:2
    for r=0:7   
        rot_ell{(k-1)*8 + r + 2} = imrotate(base_ellipse{k+1},r*180/8,'nearest','crop');
    end
end

%% ELLIPSE DATA

img = im2double(cortex_bw + cables_edge_mask);

ell_grid = zeros(size(cortex_bw,1),size(cortex_bw,2),3);

i=n+1;
while i+N-1 < size(cortex_bw,1)
   j=n+1;
   while j+N-1 < size(cortex_bw,2)
       
    block = img(i:i+N-1,j:j+N-1);
    
    cortex_count = sum(sum(double(cortex_bw(i:i+N-1,j:j+N-1))));
    cable_count = sum(sum(double(cables_edge_mask(i:i+N-1,j:j+N-1))));
    wm_count = sum(sum(double(wm(i:i+N-1,j:j+N-1))));
    
    if cortex_count > N^2/10 || (cable_count == 0 && cortex_count > 0)
        ell_grid(i:i+N-1,j:j+N-1,1:3) = rot_ell{1};
    elseif sum(sum(double(block))) > 0
        
        [~,gdir] = imgradient(block);
        block_hist = hist( reshape(gdir(gdir~=0),1,numel(gdir(gdir~=0))) , linspace(-180,180,10) );
        block_hist = block_hist + 10^-6;
        block_hist = block_hist / sum(block_hist);

        kldiv = 10^6;
        for r = 1:2*8+1
            kldiv_new = KLdiv(ell_hist{r},block_hist);
            if kldiv_new < kldiv
                kldiv = kldiv_new;
                if r== 1 && cortex_count < N^2/10
                ell_grid(i:i+N-1,j:j+N-1,1:3) = base_ellipse{4}; 
                else
                ell_grid(i:i+N-1,j:j+N-1,1:3) = rot_ell{r};
                end
            end
        end
        
    end
    j = j + N + 2*n;   
   end
    i = i + N + 2*n;
end

%% TENSOR GRID

cables_edge = color_mask(cables_edge_mask,241/255,220/255,219/255); % CABLE COLOR LIGHT 230/255,200/255,200/255
brain = color_mask(cortex_mask,218/255,150/255,149/255); % BRAIN COLOR LIGHT 210/255,170/255,170/255

img_tot = 1*brain + cables_edge + 1*cables_filling + background_color + grid;
img_tot(ell_grid~=0) = ell_grid(ell_grid~=0);
imwrite(img_tot,'output/8_ellgrid.png')


%% INSET ANISO
% inset0 = imread('inset_aniso_4.png');
% 
% inset0 = imopen(inset0,strel('disk',2));
% 
% inset0 = imresize(inset0,[600 600]); % not necessary
% 
% imwrite(inset0,'inset_aniso_5.png')
% 
% ell_inset = imread('tensors/base3.png');
% ell_inset = imresize(ell_inset,[500 500]); % 600x600 for iso
% ell_inset = imopen(ell_inset,strel('disk',10));
% ell_inset = padarray(ell_inset,[50 50]); % comment out for iso
% 
% inset1 = inset0;
% inset1(ell_inset>0) = ell_inset(ell_inset>0);
% 
% imwrite(inset1,'inset_ansio_6.png')
% 
% for i = 1:size(inset0,1)
%     for j = 1:size(inset0,2)
%         if inset0(i,j,1) > 140 && inset0(i,j,2) < 200 && inset0(i,j,3) < 150 || (inset0(i,j,1) > 200 && inset0(i,j,2) < 200)
%         inset0(i,j,1:3) = [230 200 200];
%         end
%     end
% end
% 
% imwrite(inset0,'inset_aniso_7.png')

%% INSET IMAGE

img_tot = 1*brain + cables_edge + 1*cables_filling + background_color + grid;

inset0 = imread('insets/inset_iso_7.png'); % 5 or 7

d=150;
img_tot(d:d+size(inset0,1)-1,size(img_tot,2)-size(inset0,2)+1-d:size(img_tot,2)-d,1:3) = im2double(inset0);

marked0 = imread('insets/inset_iso_2.png'); % 0 or 2

p=6;
marked0 = padarray(marked0,[p p]);
coo = [7 22]; % (6,20) for aniso, (7,22) for iso
img_tot(coo(1)*(N+2*n)-p: coo(1)*(N+2*n) + size(marked0,1)-1-p,coo(2)*(N+2*n)-p: coo(2)*(N+2*n) + size(marked0,1)-1-p,1:3) = im2double(marked0);

imwrite(img_tot,'output/5_inset.png') %enumerate different insets

%% MODEL GRID
brain_background = ones(size(brain_bw));

brain_background(wm==1) = 0;

brain_background = color_mask(brain_background,250/255,230/255,230/255);

img_tot = brain_background + grid + ell_grid + repmat(wm,[1 1 3]);
img_tot(ell_grid~=0) = ell_grid(ell_grid~=0);
imwrite(img_tot,'output/9_modelgrid.png')

























