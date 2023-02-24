function out = datanorm(X)
% Data normalization between 0 and 1 
out = X - min(X(:));
out = out ./ max(X(:))-min(X(:));
end