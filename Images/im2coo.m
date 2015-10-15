function coo = im2coo(img)

coo = zeros(numel(img),3);

m=1;
for i = 1:size(img,1)
    for j = 1:size(img,2)
        if img(i,j)>0
            coo(m,1:3) = [j, size(img,2)-i+1, 1];
            m=m+1;
        else
            coo(m,1:3) = [j, size(img,2)-i+1, 0];
            m=m+1;
        end
    end
end

S = mean(coo(:,1:2));

coo(:,1:2) = coo(:,1:2) - repmat(S,size(coo,1),1);

end