function d = distance_cramariucPLGD40(p,q,threshold)
 %threshold = -85 + 101; %% Esto es en el caso de las bases de datos de TUT 
 p1= sum( ( p - threshold ) .* (p >= threshold) .* (q == 0) );
 p2= sum( ( q - threshold ) .* (q >= threshold) .* (p == 0) );
 d = distance_cramariucLGD(p,q) + (1/40)*(p1+p2);
 return
end