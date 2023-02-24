function distances = distance_neyman(p,q,~)
if max([p(:),q(:)]) <= 1.0001
factorzero = 0.0001/104;
else
factorzero = 0.0001;
end

divisor   = p + (factorzero*(p==0));
distances = sum((((p-q).^2)./divisor));

return