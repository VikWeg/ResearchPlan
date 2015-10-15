function img = coo2im(coo,Nx,Ny)

maxX = max(coo(:,1));
maxY = max(coo(:,2));

img = zeros(Nx,Ny);

for i = 1:size(coo,1)
    img( Ny - round(Ny*(coo(i,2)+maxY+1)/(2*maxY+1)) + 2, round((Nx-2)*(coo(i,1)+maxX+1)/(2*maxX+1)) + 1 ) = 1;
end

img = padarray(img,[10 10]);

img = imdilate(img,strel('disk',3));

img_edge = edge(img,'Canny');
img_edge_dilated = imdilate(img_edge,strel('disk',1));
img_filtered = imgaussfilt(img,1);
img(img_edge_dilated) = img_filtered(img_edge_dilated);

end