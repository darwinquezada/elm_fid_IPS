close all; clear all; clc;

load('/Users/darwinquezada/Documents/MATLAB/DataCompression/databases/LIB1C.mat');
addpath('databases','elm_functions','read_functions','positioning_functions','Results','clusters');
%Initialization
k=1;
Eps=7;
MinPts=5;
datarep='positive';
distanceMetric='distance_euclidean';
labelk2=1;
databaseName='LIB1C';
clustering_method = 'knn_dbscan';

[IDX] = dbscan(database.trainingMacs,Eps,MinPts,'Distance','euclidean');

sep = filesep();

%Computing DBSCAN clustering
fprintf('Starting clustering DBSCAN');
    if strcmp(clustering_method,'knn_dbscan')   
        % knn_dbscan( database_orig , k, datarep, defaultNonDetectedValue,
        % newNonDetectedValue ,distanceMetric , k2 , randseed, labelk2, folderKMeansBase)
        
        folderClusters = ['clusters' sep databaseName sep 'knn_dbscan' sep] ;
        
        [results]  = feval(clustering_method, database , k, Eps, MinPts, datarep, ...
            distanceMetric, 1, labelk2, folderClusters); 
        distance = strsplit(distanceMetric,'_');
        folderResults = ['Results' sep databaseName sep clustering_method sep  datarep '_' distance{2} '_k' sprintf('%03d',k) sep] ;
    end
    
    disp(folderResults);
    
    if ~exist(folderResults)
      mkdir(folderResults)
    end
    save([folderResults 'results_rep.mat' ],'results'); % Modificado por Darwin
