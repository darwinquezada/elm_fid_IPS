function [x,y,z,error2D,error3D] = error_estimation(estimatedlabels, testLabels)
    rowEstimatedLabels = size(estimatedlabels,1);
    dist_3D = zeros(rowEstimatedLabels,1);
    dist_2D = zeros(rowEstimatedLabels,1);
    for i = 1:rowEstimatedLabels
        dist_3D(i) = sqrt(sum((estimatedlabels(i,1:3)-testLabels(i,1:3)).^2));
        dist_2D(i) = sqrt(sum((estimatedlabels(i,1:2)-testLabels(i,1:2)).^2));
    end
    x = 0;
    y = 0;
    z = 0;
    error2D = mean(dist_2D);
    error3D = mean(dist_3D);
    
end