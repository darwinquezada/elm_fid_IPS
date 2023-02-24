function distances = distance_g3_sorensen(p,q,~)

distances = sum(abs(p-q)) / sum(p+q);

return