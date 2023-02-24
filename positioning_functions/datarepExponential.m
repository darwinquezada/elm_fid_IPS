function [ db1 ] = datarepExponential( db0 )

    minValue           = min([db0.trainingMacs(:)',db0.testMacs(:)']);

    normValue          = exp(-minValue/24);
    
    db1.trainingMacs   = exp((db0.trainingMacs-minValue)/24)./normValue;
    db1.testMacs       = exp((db0.testMacs-minValue)/24)./normValue;

    db1.trainingLabels   = db0.trainingLabels;    
    db1.testLabels       = db0.testLabels;
    
end
