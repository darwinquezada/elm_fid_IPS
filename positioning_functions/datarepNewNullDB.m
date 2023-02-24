function [ db1 ] = datarepNewNullDB( db0 , oldNull, newNull )

    db1.trainingMacs   = datarepNewNull( db0.trainingMacs  , oldNull, newNull );
    db1.testMacs       = datarepNewNull( db0.testMacs      , oldNull, newNull );

    db1.trainingLabels   = db0.trainingLabels;    
    db1.testLabels       = db0.testLabels;
    
end
