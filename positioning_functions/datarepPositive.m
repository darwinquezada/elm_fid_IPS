function [ db1 ] = datarepPositive( db0 )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

    minValue           = min([db0.trainingMacs(:)',db0.testMacs(:)']);
    
    db1.trainingMacs   = db0.trainingMacs-minValue;
    db1.testMacs       = db0.testMacs-minValue;

    db1.trainingLabels = db0.trainingLabels;
    db1.testLabels     = db0.testLabels;

end
