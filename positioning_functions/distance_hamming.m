function distances = distance_hamming(p,q,~)

distances = sum(p~=q);

return