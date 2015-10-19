function masked_rgb = mask_rgb(im, mask)

r = im(:,:,1);
g = im(:,:,2);
b = im(:,:,3);

r = r.*mask;
g = g.*mask;
b = b.*mask;

masked_rgb = zeros(size(im));

masked_rgb(:,:,1) = r;
masked_rgb(:,:,2) = g;
masked_rgb(:,:,3) = b;

end