function masked_rgb = color_mask(mask,rval,gval,bval)

dim = size(mask);

masked_rgb = zeros(dim(1),dim(2),3);

masked_rgb(:,:,1) = rval*mask;
masked_rgb(:,:,2) = gval*mask;
masked_rgb(:,:,3) = bval*mask;

end