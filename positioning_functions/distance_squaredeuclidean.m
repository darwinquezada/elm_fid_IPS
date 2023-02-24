function distances = distance_squaredeuclidean(p,q,~)

distances = sum((p-q).^2);

return