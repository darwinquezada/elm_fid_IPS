function [ results , clusters ] = knn_dbscan( database_orig , k, Eps, MinPts, datarep, distanceMetric, randseed, labelk2, folderDbscanBase)
% knn_ips Indoor Positioning System based on kNN
%   
% inputs: 
%    database_orig           : database with RSSI values and Labels
%    k                       : value of nearest neighbors in kNN algoritm
%    datarep                 : Data representation
%    defaultNonDetectedValue : Value in DB for non detected RSSI values
%    newNonDetectedValue     : New Value for non detected RSSI values
%    distanceMetric          : Distance metric used for FP matching
%
% output:
%    results                 : output statistics
%
% Developed by J. Torres-Sospedra,
% Instiute of New Imaging Technologies, Universitat Jaume I
% jtorres@uji.es
% Modified by Darwin Quezada

defaultNonDetectedValue = 100; 
newNonDetectedValue     = min([database_orig.trainingMacs(:);database_orig.testMacs(:)]) - 1;

ticT1m3 = myTic();
% Remap bld numbers to 1:nblds
origBlds = unique((database_orig.trainingLabels(:,5))); 
nblds    = size(unique(origBlds),1);
database0  = remapBldDB(database_orig,origBlds',1:nblds);

% Remap floors numbers to 1:nfloors
origFloors = unique((database_orig.trainingLabels(:,4))); 
nfloors    = size(unique(origFloors),1);
database0  = remapFloorDB(database0,origFloors',1:nfloors);

% Calculate the overall min RSSI (minus 1)
if size(newNonDetectedValue,1) == 0
newNonDetectedValue = min(...
                     [database0.trainingMacs(:)',...                      
                      database0.testMacs(:)']...
                     )-1;
end

% Replace NonDetectedValue
if size(defaultNonDetectedValue,1) ~= 0
    database0 = datarepNewNullDB(database0,defaultNonDetectedValue,newNonDetectedValue);
end

% Apply new data representation                 
if  strcmp(datarep,'positive') % Convert negative values to positive by adding a value
    database  = datarepPositive(database0);
  
    if or(strcmp(distanceMetric,'distance_cramariucPLGD40'),strcmp(distanceMetric,'distance_cramariucPLGD40'))
      additionalparams = -85 - newNonDetectedValue;
    else
      additionalparams = 0;
    end
    
elseif strcmp(datarep,'exponential') % Convert to exponential data representation
    %database0 = datarepPositive(database0);
    database  = datarepExponential(database0);
    
    if or(strcmp(distanceMetric,'distance_cramariucPLGD40'),strcmp(distanceMetric,'distance_cramariucPLGD40'))
      additionalparams =  exp(((-85-newNonDetectedValue))/24)./exp((-newNonDetectedValue)/24);
    else
      additionalparams = 0;
    end
    
elseif strcmp(datarep,'powed') % Convert to powed data representation
    %database0 = datarepPositive(database0); 
    database  = datarepPowed(database0);  
    
    if or(strcmp(distanceMetric,'distance_cramariucPLGD40'),strcmp(distanceMetric,'distance_cramariucPLGD40'))
      additionalparams =  ((-85-newNonDetectedValue).^exp(1))./((-newNonDetectedValue).^exp(1));
    else
      additionalparams = 0;
    end
      
end % Default, no conversion

clear database0;
   
% Get number of samples
rsamples = size(database.trainingMacs,1); 
osamples = size(database.testMacs,1);
nmacs    = size(database.testMacs,2);

% Prepare output vectors
results.error      = zeros(osamples,5);
results.prediction = zeros(osamples,5);
results.targets    = zeros(osamples,5);
results.candidates = zeros(osamples,1);
results.distances  = zeros(osamples,1);
results.timesample = zeros(osamples,5);
results.considered = zeros(osamples,5);
results.cluster    = zeros(osamples,1);
results.csamples   = zeros(osamples,1);

fprintf('\n')
fprintf('    database features      : [%d,%d,%d]\n',rsamples,osamples,nmacs)
fprintf('    k                      : %d\n',k)
fprintf('    Eps                    : %d\n',Eps)
fprintf('    MinPts                 : %d\n',MinPts)
fprintf('    datarep                : %s\n',datarep)
fprintf('    defaultNonDetectedValue: %d\n',defaultNonDetectedValue)
fprintf('    newNonDetectedValue    : %d\n',newNonDetectedValue)
fprintf('    distanceMetric         : %s\n',distanceMetric )

time00pre = myToc(ticT1m3);

database.trainingValidMacs = (database_orig.trainingMacs ~= defaultNonDetectedValue); 
database.testValidMacs     = (database_orig.testMacs     ~= defaultNonDetectedValue);
% vecidxmacs    = 1:nmacs;
vecidxsamples = 1:rsamples;

sprintf('%03d',randseed);

file = [folderDbscanBase,'clusters_dbscan_',datarep,'_',...
    sprintf('%03d',randseed),'.mat'];

fileClusters = strcat(file);

if ~exist(fileClusters)
	clusters.timeDbscan01 = myToc(ticT1m3);
	if ~exist(folderDbscanBase)
		mkdir(folderDbscanBase);
	end
	clusters.randseed = rsamples+osamples+randseed;
	rand('seed',clusters.randseed);

    distance = strsplit(distanceMetric,'_');
    
    % Init DBSCAN
    distance = strsplit(distanceMetric,'_');
    [clusters.idx,clusters.noise,clusters.centroids,clusters.sumd] = local_dbscan(database.trainingMacs,Eps,MinPts,'cityblock'); 

    % End DBSCAN
    
	clusters.timeDbscan02 = myToc(ticT1m3);
	
    save(fileClusters,'clusters');
else 
    load(fileClusters,'clusters'); 
end

csamples = size(clusters.centroids,1); 
idxs     = cell(1,csamples);
for i = 1:csamples
    idxs{i} = vecidxsamples(clusters.idx == i);
end

time00 = myToc(ticT1m3);

srt = sort(clusters.idx);
GC = zeros(max(srt)+1,1);
GR = zeros(max(srt)+1,1);
b=1;
for i = 0:max(srt)
    GR(b) = i;
    for x = 1:size(srt,1)
        if GR(b) == srt(x)
            GC(b) = GC(b) + 1;
        end
    end
    b=b+1;
end

numgr = size(GR(GR~=0,:));
fprintf('\nClustering \n')
fprintf('    Eph for dbscan             : %d\n',Eps)
fprintf('    Minpts for dbscan          : %d\n',MinPts)
fprintf('    number of clusters         : %d\n',numgr(1))
fprintf('    Noise (Denoted by 0 or -1) : %d\n\n',GC(GR==0,:))

for i = 1:osamples  

    time01 = myToc(ticT1m3);
        
    ofp = database.testMacs(i,:);
	
    dist_cluster = zeros(1,csamples);
		%% Find Cluster
    for j = 1:csamples
      dist_cluster(1,j) = distance_cityblock(clusters.centroids(j,:),ofp);
    end
    [value,cluster] = min(dist_cluster);
    
	samples = idxs{cluster};
	
    if size(samples,1)==0
      samples = 1:rsamples;      
    end
  
	  trainingMacs      = database.trainingMacs(samples,:);
	  trainingValidMacs = database.trainingValidMacs(samples,:);	
	  trainingLabels    = database.trainingLabels(samples,:);
	
	  rsamplesreduced = size(trainingMacs,1);
    
    ofp = database.testMacs(i,:);
    
    time02 = myToc(ticT1m3);
    
    distances = zeros(1,rsamplesreduced);
 
    for j = 1:rsamplesreduced
        distances(j) = feval(distanceMetric,trainingMacs(j,:),ofp,additionalparams);
    end
           
    time03 = myToc(ticT1m3);
    % Sort distances
    [distancessort,candidates] = sort(distances);
    
    % Set the value of candidates to k (kNN)  <---- THIS IS NEW!!!!!!
    % MENTION IN THE PAPER
    ncandidates = min(k,rsamplesreduced);
    
    % Check if candidates ranked in positions higher than k have the same
    % distance in the feature space.
    while ncandidates < rsamplesreduced 
       if abs(distancessort(ncandidates)-distancessort(ncandidates+1))<0.000000000001
           ncandidates = ncandidates+1;
       else
           break
       end
    end

    % Estimate the building from the ncandidates
    probBld = zeros(1,nblds);
    for bld = 1:nblds
        probBld(bld) = sum(trainingLabels(candidates(1:ncandidates),5)==bld);
    end
    [~,bld] = max(probBld);
    
    % Estimate the floor from the ncandidates
    probFloor = zeros(1,nfloors);
    for floor = 1:nfloors
        probFloor(floor) = sum( (trainingLabels(candidates(1:ncandidates),4)==floor).*...
                                (trainingLabels(candidates(1:ncandidates),5)==bld));
    end
    [~,floor] = max(probFloor);
    
    % Estimate the coordinates from the ncandidates that belong to the
    % estimated floor
    points = 0;
    point  = [0,0,0];
    for j = 1:ncandidates
        if (trainingLabels(candidates(j),4)==floor).*(trainingLabels(candidates(j),5)==bld)
            points = points + 1;
            point  = point + trainingLabels(candidates(j),1:3);
        end
    end
    point = point / points;
    
    time04 = myToc(ticT1m3);
    % Extract the real-world bld (revert bld remap)
    realWorldBld = remapVector(bld,1:nblds,origBlds');
    
    % Extract the real-world floor (revert floor remap)
    realWorldFloor = remapVector(floor,1:nfloors,origFloors');
    
    % Generate the results statistics
    results.error(i,1) = distance_euclidean(database.testLabels(i,1:2),point(1:2)); % 2D Positioning error in m 
    results.error(i,2) = distance_euclidean(database.testLabels(i,1:3),point(1:3)); % 3D Positioning error in m
    results.error(i,3) = abs(database.testLabels(i,3)-point(3));                    % Height difference in m
    results.error(i,4) = abs(realWorldFloor - database_orig.testLabels(i,4));       % Height difference in floors
    results.error(i,5) = (realWorldBld ~= database_orig.testLabels(i,5));           % Building error
    
    results.prediction(i,:) = [point,realWorldFloor,realWorldBld];                  % Predicted position [x,y,z,floor]
    results.targets(i,:)    = database_orig.testLabels(i,1:5);                      % Current position   [x,y,z,floor]
    results.candidates(i,1) = ncandidates;                                          % Number of nearest neighbours used to estimate the position
    results.distances(i,1)  = distancessort(1);                                     % Distance in feature space of the best match
    
	results.cluster(i,1)    = cluster;
	results.csamples(i,1)   = csamples;	
	
    results.samples(i,1)    = rsamplesreduced;         % Distance in feature space of the best match
    results.timesample(i,1) = time01;                  % Distance in feature space of the best match
    results.timesample(i,2) = time02;                  % Distance in feature space of the best match
    results.timesample(i,3) = time03;                  % Distance in feature space of the best match    
    results.timesample(i,4) = time04;                  % Distance in feature space of the best match    
    results.timesample(i,5) = myToc(ticT1m3);                     % Distance in feature space of the best match    
    
    time05 = myToc(ticT1m3);
	
end

results.timefull = [time00pre,time00,myToc(ticT1m3)];

fprintf('Successfully executed!\n')
fprintf('    Mean Positioning Error 2D*  : %3.3f\n',mean(results.error((results.error(:,4)==0),1))*sqrt(0.83673469));  % Mean 2D Positioning error when estimate and true position are in the same floor
fprintf('    Mean Positioning Error 3D   : %3.3f\n',mean(results.error(:,2))*sqrt(0.83673469));   
% fprintf('    Mean Positioning Error EvAAL: %3.3f\n',mean(results.error(:,1)+4*results.error(:,4)+50*results.error(:,5))*sqrt(0.83673469));                        % Mean 3D Positioning error
% fprintf('    Mean Positioning Error IPIN : %3.3f\n',mean(results.error(:,1)+15*results.error(:,4)+50*results.error(:,5))*sqrt(0.83673469));                        % Mean 3D Positioning error
fprintf('    Bld&Floor detection hit rate: %3.3f\n',mean((results.error(:,4)==0).*(results.error(:,5)==0))*100);                 % Floor detection rate
fprintf('    BD Preprocessing time       : %3.6f\n',time00);    
fprintf('    Average preprocessing time  : %3.6f\n',mean(results.timesample(:,2)-results.timesample(:,1)));    
fprintf('    Average matching time       : %3.6f\n',mean(results.timesample(:,3)-results.timesample(:,2)));    
fprintf('    Average positioning time    : %3.6f\n',mean(results.timesample(:,4)-results.timesample(:,3)));   
fprintf('    Full time                   : %3.6f\n',results.timefull(3)); 
fprintf('        DB Preprocessing        : %3.3f%%\n',results.timefull(1)/results.timefull(3)*100); 
fprintf('        FP Preprocessing        : %3.3f%%\n',sum(results.timesample(:,2)-results.timesample(:,1))/results.timefull(3)*100); 
fprintf('        Matching                : %3.3f%%\n',sum(results.timesample(:,3)-results.timesample(:,2))/results.timefull(3)*100); 
fprintf('        Positioning             : %3.3f%%\n\n',sum(results.timesample(:,4)-results.timesample(:,3))/results.timefull(3)*100); 