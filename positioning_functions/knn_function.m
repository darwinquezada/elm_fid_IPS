function estimatedlabels = knn_function(dataTraining,dataTest,dataTrainingLabels,distance,k)
    
    distance = distance(find(~isspace(distance)));
    dmatrix = distance_matrix(dataTraining, dataTest, distance);
    rowTest = size(dataTest,1);
    estimatedlabels = zeros(1,3);
    for i = 1:rowTest
        vpos = zeros(k,1); 
        for j = 1:k
            [vval vpos(j)] = min(dmatrix(i,:));
            dmatrix(i,vpos(j)) = 1e6;
        end

        x = 0; y = 0; z = 0;

        for l = 1:k
          x  = x  + dataTrainingLabels(vpos(l),1);
          y = y + dataTrainingLabels(vpos(l),2);
          z = z + dataTrainingLabels(vpos(l),3);
        end

        estimatedlabels(i,:) = [x/k, y/k, z/k];    
    end
end