function W = initializeSVDO(neurons,trainingMacs)
    [U,S,V] = svd(trainingMacs);
    W = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))).';
    % W = normalize(V);
    W = W(1:neurons,:);