function [ db1 ] = remapBldDB( db0 , bldin, bldout )

    db1.trainingMacs   = db0.trainingMacs;
    db1.testMacs       = db0.testMacs;

    db1.trainingLabels   = db0.trainingLabels;
    db1.testLabels       = db0.testLabels;
        
    for i = 1:size(bldin,2)
        db1.trainingLabels(db0.trainingLabels(:,5)==bldin(i),5) = bldout(i);
        db1.testLabels(db0.testLabels(:,5)==bldin(i),5)         = bldout(i);
    end
    
end
