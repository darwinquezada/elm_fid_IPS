function [ db1 ] = datarepPowed( db0 )
  
    minValue           = min([db0.trainingMacs(:)',db0.testMacs(:)']);
    
    normValue          = ((-minValue).^exp(1));
    
    db1.trainingMacs   = ((db0.trainingMacs-minValue).^exp(1))./normValue;
    db1.testMacs       = ((db0.testMacs-minValue).^exp(1))./normValue;

    db1.trainingLabels   = db0.trainingLabels;
    db1.testLabels       = db0.testLabels;
    
end
