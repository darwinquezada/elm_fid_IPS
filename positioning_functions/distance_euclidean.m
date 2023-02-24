function distances = distance_euclidean(p,q,~)

distances = sqrt(sum((p-q).^2));

return