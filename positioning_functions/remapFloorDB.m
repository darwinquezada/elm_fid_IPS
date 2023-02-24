function [ db1 ] = remapFloorDB( db0 , floorin, floorout )

    db1.trainingMacs   = db0.trainingMacs;
    db1.testMacs       = db0.testMacs;

    db1.trainingLabels   = db0.trainingLabels;
    db1.testLabels       = db0.testLabels;
        
    for i = 1:size(floorin,2)
        db1.trainingLabels(db0.trainingLabels(:,4)==floorin(i),4)     = floorout(i);
        db1.testLabels(db0.testLabels(:,4)==floorin(i),4)             = floorout(i);
    end
    
end
