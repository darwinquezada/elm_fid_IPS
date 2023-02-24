function distances = distance_cityblock(p,q,~)

distances = sum(abs(p-q));

return