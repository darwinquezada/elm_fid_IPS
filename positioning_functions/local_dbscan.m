%
% Copyright (c) 2015, Yarpiz (www.yarpiz.com)
% All rights reserved. Please read the "license.txt" for license terms.
%
% Project Code: YPML110
% Project Title: Implementation of DBSCAN Clustering in MATLAB
% Publisher: Yarpiz (www.yarpiz.com)
% 
% Developer: S. Mostapha Kalami Heris (Member of Yarpiz Team)
% 
% Contact Info: sm.kalami@gmail.com, info@yarpiz.com
%
% Last modification : DQ
%       Distance :
%         'euclidean'        - Euclidean distance (default)
%         'squaredeuclidean' - Squared Euclidean distance
%         'seuclidean'       - Standardized Euclidean distance. Each
%                              coordinate difference between rows in X and Y
%                              is scaled by dividing by the corresponding
%                              element of the standard deviation computed
%                              from X, S=NANSTD(X). To specify another value
%                              for S, use
%                              D = pdist2(X,Y,'seuclidean',S).
%         'cityblock'        - City Block distance
%         'minkowski'        - Minkowski distance. The default exponent is 2.
%                              To specify a different exponent, use
%                              D = pdist2(X,Y,'minkowski',P), where the
%                              exponent P is a scalar positive value.
%         'chebychev'        - Chebychev distance (maximum coordinate
%                              difference)
%         'mahalanobis'      - Mahalanobis distance, using the sample
%                              covariance of X as computed by NANCOV.  To
%                              compute the distance with a different
%                              covariance, use
%                              D = pdist2(X,Y,'mahalanobis',C), where the
%                              matrix C is symmetric and positive definite.
%         'cosine'           - One minus the cosine of the included angle
%                              between observations (treated as vectors)
%         'correlation'      - One minus the sample linear correlation
%                              between observations (treated as sequences of
%                              values).
%         'spearman'         - One minus the sample Spearman's rank
%                              correlation between observations (treated as
%                              sequences of values)
%         'hamming'          - Hamming distance, percentage of coordinates
%                              that differ
%         'jaccard'          - One minus the Jaccard coefficient, the
%                              percentage of nonzero coordinates that differ
function [IDX, isnoise, centroid, sumd]=local_dbscan(X,epsilon,MinPts,Distance)
    C=0;
    
    n=size(X,1);
    IDX=zeros(n,1);
    
    D=distance_matrix(X,X,Distance);  
    
    visited=false(n,1);
    isnoise=false(n,1);
    
    for i=1:n
        if ~visited(i)
            visited(i)=true;
            
            Neighbors=RegionQuery(i);
            if numel(Neighbors)<MinPts
                % X(i,:) is NOISE
                isnoise(i)=true;
            else
                C=C+1;
                ExpandCluster(i,Neighbors,C);
            end
            
        end
    
    end
    
    function ExpandCluster(i,Neighbors,C)
        IDX(i)=C;
        
        k = 1;
        while true
            j = Neighbors(k);
            
            if ~visited(j)
                visited(j)=true;
                Neighbors2=RegionQuery(j);
                if numel(Neighbors2)>=MinPts
                    Neighbors=[Neighbors Neighbors2];   %#ok
                end
            end
            if IDX(j)==0
                IDX(j)=C;
            end
            
            k = k + 1;
            if k > numel(Neighbors)
                break;
            end
        end
    end
    
    function Neighbors=RegionQuery(i)
        Neighbors=find(D(i,:)<=epsilon);
    end

    srt = sort(IDX);
    SC = zeros(max(srt)+1,1);
    C = zeros(max(srt)+1,1);
    b=1;
    for i = 0:max(srt)
        C(b) = i;
        for x = 1:size(srt,1)
            if C(b) == srt(x)
                SC(b) = SC(b) + 1;
            end
        end
        b=b+1;
    end
    % [SC,C] = groupcounts(IDX); %Group the results
    % SC = number of samples per cluster
    % C = Clusters
    rws = size(X,1);
    cols = size(X,2);
    max_num_clusters = max(C);

    sc_no_noise = SC(1:size(SC),1);
    sumd = zeros(max_num_clusters,1);
    % Look for the centroids
    k = 1; 
    centroid = zeros(max_num_clusters,cols);
    distc = zeros(max_num_clusters,cols);
    for i = 1:max_num_clusters
        inc = 1;
        new_cluster = zeros(sc_no_noise(i),cols);
        for j = 1:rws
            if IDX(j)==i
                new_cluster(inc,:)= X(j,:);
                inc = inc+1;
            end
        end
        [~,centroid(i,:),sumd(i)] = kmeans(new_cluster,k);
    end
end