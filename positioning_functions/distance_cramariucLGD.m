function d = distance_cramariucLGD(p,q)
 sigma       = 5;
 threshold   = 0.0001;
 
 numerator   = -((p - q).^2);
 denominator = 2*(sigma^2);
 differences = (1/sqrt(2*pi*sigma^2)).*exp(numerator/denominator);
 nonzero     = (p.*q)~=0;
 d = -sum( log (differences.*nonzero + threshold*(1-nonzero)) );
 
end