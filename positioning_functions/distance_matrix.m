function disMatrix = distance_matrix(dataTraining, dataTest, distance)
    
    rowTest = size(dataTest,1);
    rowTraining = size(dataTraining,1);
    disMatrix = zeros(rowTest,rowTraining);
    distance = distance(find(~isspace(distance)));

    for i = 1:rowTest
        for j = 1:rowTraining
            switch lower(distance)
                case 'euclidean'
                    disMatrix(i,j) = distance_euclidean(dataTest(i,:),dataTraining(j,:));
                case 'squaredeuclidean'
                    disMatrix(i,j) = distance_squaredeuclidean(dataTest(i,:),dataTraining(j,:));
                case 'manhattan'
                    disMatrix(i,j) = manhattan(dataTest(i,:),dataTraining(j,:));
                case 'minkowsky'
                    disMatrix(i,j) = minkowsky(dataTest(i,:),dataTraining(j,:));
                case 'chevyshev'
                    disMatrix(i,j) = chevyshev(dataTest(i,:),dataTraining(j,:));
                case 'neyman'
                    disMatrix(i,j) = neyman(dataTest(i,:),dataTraining(j,:));
                case 'cityblock'
                    disMatrix(i,j) = distance_cityblock(dataTest(i,:),dataTraining(j,:));
                case 'hamming'
                    disMatrix(i,j) = distance_hamming(dataTest(i,:),dataTraining(j,:));
                case 'g3_sorensen'
                    disMatrix(i,j) = distance_g3_sorensen(dataTest(i,:),dataTraining(j,:));
                case 'cramariucplgd10'
                    disMatrix(i,j) = distance_cramariucPLGD10(dataTest(i,:),dataTraining(j,:));
                otherwise
                  disp('Unknown distance.')
            end
        end
    end
end

%% Manhattan distance
function distance = manhattan(vtest, vtraining)
    distance = sum(abs(vtest-vtraining));
end

%% Minkowsky
function distance = minkowsky(vtest, vtraining)
    distance = (sum(abs(vtest-vtraining).^5))^(1/5);
end

%% Chevyshev
function distance = chevyshev(vtest, vtraining)
    distance = max(abs(vtest-vtraining));
end

%% Neyman
function distance = neyman(vtest, vtraining)
  if max([vtest(:),vtraining(:)]) <= 1.0001
    factor = 0.0001/104;
  else
    factor = 0.0001;
  end
  divisor   = vtraining + (factor*(vtraining==0));
  distance = sum((((vtest-vtraining).^2)./divisor));
end