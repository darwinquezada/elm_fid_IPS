function database = get_data_representation(database_orig, datarep, distanceMetric)
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
elseif strcmp(datarep,'posuint') % Convert to powed data representation
    %database0 = datarepPositive(database0); 
    database  = datarepPosuint(database0);  
    
    if or(strcmp(distanceMetric,'distance_cramariucPLGD40'),strcmp(distanceMetric,'distance_cramariucPLGD40'))
      additionalparams =  ((-85-newNonDetectedValue).^exp(1))./((-newNonDetectedValue).^exp(1));
    else
      additionalparams = 0;
    end      
end % Default, no conversion
