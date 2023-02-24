function [mean2d, mean3d, fhr, bhr] = svm_position(database)
    % Floor hit rate
    model_floor=fitcsvm(database.trainingMacs, database.trainingLabels(:,4));
    predict_floor=svm.predict(model_floor,database.testMacs);
    display(predict_floor)
    % Building hit rate
    addpath("positioning_functions")

end 