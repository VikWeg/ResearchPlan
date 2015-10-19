function ell = im2ell(img)

nx = size(img,1);
ny = size(img,2);

coo = im2coo(img);

[~,S,V] = svd(coo(logical(coo(:,3)),1:2));

t = linspace(0,2*pi);

ell_coo = [sqrt(S(1,1))*cos(t)', sqrt(S(2,2))*sin(t + asin(V(2,1)))'];

ell = coo2im(ell_coo,nx,ny);

ell = ell(1:nx,1:ny);

end