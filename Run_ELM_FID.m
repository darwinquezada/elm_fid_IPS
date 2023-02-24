% Running all experiments 

clear all; close all; clc;

addpath('databases','elm_functions','read_function','positioning_functions', 'Initializers');


% Databases
databases = {'DSI1','DSI2','LIB1','LIB2', 'MAN1', 'MAN2', 'TUT1','TUT2','TUT3','TUT4','TUT5','TUT6','TUT7'};
% Final database Dimension
dimensions = [46,43,31,45,15,8,50,55,113,182,6,53,116];
% Data representation
representation = {'positive','positive','positive','positive','positive','positive','positive','positive','positive',...
    'positive','positive','positive','positive'};
% k value (kNN)
k = [1,1,1,1,1,1,1,1,1,1,1,1,1];
% Distance metric
distance = {'cityblock','cityblock','cityblock','cityblock','cityblock','cityblock','cityblock',...
    'cityblock','cityblock','cityblock','cityblock','cityblock','cityblock'};
% Activation functions
act_functions = { 'sin','cos','sig','tansig','relu' };

defaultNonDetectedValue = 100; 
mu = 0; % mean 
sigma = 0.01; % Standard deviation

report = zeros(length(databases)* length(act_functions),8);

% Number of iterations
iterations = 20;
row_count = 0;

for idx = 1:length(databases)
    for i = 1:length(act_functions)
        row_count = row_count + 1;
        disp(strcat("--------- Database ", databases{idx}," ----------------"))
        % Load dataset
        load(strcat('databases/',databases{idx},'.mat'));
        newNonDetectedValue = min([database.trainingMacs(:);database.testMacs(:)]) - 1;
        if size(defaultNonDetectedValue,1) ~= 0
            database = datarepNewNullDB(database,defaultNonDetectedValue,newNonDetectedValue);
        end

        % Change data representation
        db0 = data_representation(database,representation{idx});
        trainingMacs = db0.trainingMacs;
        testMacs = db0.testMacs;

        trainingLabels = db0.trainingLabels;
        testLabels = db0.testLabels;

        neurons = dimensions(idx);

        % FID
        [U,S,V] = svd(trainingMacs);
        AF = act_functions{i};                   % Activation Function
        W = bsxfun(@rdivide, V, sqrt(sum(V.^2, 2))).';
        W = W(1:neurons,:);

        pos_error = zeros(iterations,7);
        floor_rate = zeros(iterations,7);
        time_initialization = zeros(iterations,7);
        time_compression_train = zeros(iterations,7);
        time_compression_test = zeros(iterations,7);
        time_positioning = zeros(iterations,7);

        features = size(trainingMacs);

        for iter = 1:iterations

            WHe = initializeHe([neurons size(trainingMacs,2)], 1);
            WGau = initializeGaussian([size(trainingMacs,2) neurons],mu,sigma);
            Wort = initializeOrthogonal([size(trainingMacs,2) neurons]);
            WGlo = initializeGlorot([size(trainingMacs,2) neurons],1,1);
            WUni = initializeUniform([size(trainingMacs,2) neurons],0.1);
            time01 = cputime;
            r = randperm(features(1), neurons); 
            Winp = normalize(trainingMacs(r,:),'zscore');
            Winp(isnan(Winp))=0;
            timeIMP = cputime-time01;
            time_initialization(iter,1)=timeIMP;

            %% Orthogonal initialization
            disp('----------------- Orthogonal Initialization ---------------------');
            [Bort_train,ZHort_train,Hort_train,RMSEort_train, Wort_train] = ae_elm_encoder(trainingMacs,AF,neurons);
            [Bort_test,ZHort_test,Hort_test,RMSEort_test] = ae_elm_encoder_w(testMacs,AF,neurons,Wort_train);
    
            databaseOrt.trainingMacs = Hort_train;
            databaseOrt.testMacs = Hort_test;
            databaseOrt.trainingLabels = database.trainingLabels;
            databaseOrt.testLabels = database.testLabels;
    
            [ results_kNNOrt ] = knn_baseline( databaseOrt , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'Orthogonal/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNOrt');

            pos_error(iter,1) = mean(results_kNNOrt.error(:,2));
            floor_rate(iter,1) = mean((results_kNNOrt.error(:,4)==0).*(results_kNNOrt.error(:,5)==0))*100;

            
            %% HE initialization
            disp('----------------- He Initialization ---------------------');
            [BHE_train,ZHHE_train,HHE_train,RMSEHE_train,WHE_train]=ae_elm_encoder_w(trainingMacs,AF,neurons,WHe);   % Autoencoder
            [BHE_test,ZHHE_test,HHE_test,RMSEHE_test]=ae_elm_encoder_w(testMacs,AF,neurons,WHE_train);
    
            databaseHe.trainingMacs = HHE_train;
            databaseHe.testMacs = HHE_test;
            databaseHe.trainingLabels = database.trainingLabels;
            databaseHe.testLabels = database.testLabels;
    
            [ results_kNNHe ] = knn_baseline( databaseHe , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'HE/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNHe');

            pos_error(iter,2) = mean(results_kNNHe.error(:,2));
            floor_rate(iter,2) = mean((results_kNNHe.error(:,4)==0).*(results_kNNHe.error(:,5)==0))*100;
    
            %% Glorot initializer
            disp('----------------- Glorot Initialization ---------------------');
            [BGlo_train,ZHGlo_train,HGlo_train,RMSEGlo_train, WGlo_train]=ae_elm_encoder_w(trainingMacs,AF,neurons,WGlo');   % Autoencoder
            [BGlo_test,ZHGlo_test,HGlo_test,RMSEGlo_test]=ae_elm_encoder_w(testMacs,AF,neurons,WGlo_train);   % Autoencoder
    
            databaseGlo.trainingMacs = HGlo_train;
            databaseGlo.testMacs = HGlo_test;
            databaseGlo.trainingLabels = database.trainingLabels;
            databaseGlo.testLabels = database.testLabels;
    
            [ results_kNNGlo ] = knn_baseline( databaseGlo , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'Glorot/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNGlo');

            pos_error(iter,3) = mean(results_kNNGlo.error(:,2));
            floor_rate(iter,3) = mean((results_kNNGlo.error(:,4)==0).*(results_kNNGlo.error(:,5)==0))*100;
    
            %% Gaussian initializer
            disp('----------------- Gaussian Initialization ---------------------');
            [Bgau_train,ZHgau_train,Hgau_train,RMSEgau_train,Wgau_train]= ae_elm_encoder_w(trainingMacs,AF,neurons,WGau');   % Autoencoder
            [Bgau_test,ZHgau_test,Hgau_test,RMSEgau_test,]= ae_elm_encoder_w(testMacs,AF,neurons,Wgau_train);
    
            databaseGau.trainingMacs = Hgau_train;
            databaseGau.testMacs = Hgau_test;
            databaseGau.trainingLabels = database.trainingLabels;
            databaseGau.testLabels = database.testLabels;
    
            [ results_kNNGau ] = knn_baseline( databaseGau , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'Gaussian/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNGau');

            pos_error(iter,4) = mean(results_kNNGau.error(:,2));
            floor_rate(iter,4) = mean((results_kNNGau.error(:,4)==0).*(results_kNNGau.error(:,5)==0))*100;
    
            %% Uniform initializer
            disp('----------------- Uniform Initialization ---------------------');
            [BUni_train,ZHUni_train,HUni_train,RMSEUni_train, WUni_train]=ae_elm_encoder_w(trainingMacs,AF,neurons,WUni');   % Autoencoder
            [BUni_test,ZHUni_test,HUni_test,RMSEUni_test]=ae_elm_encoder_w(testMacs,AF,neurons,WUni_train);   % Autoencoder
    
            databaseUni.trainingMacs = HUni_train;
            databaseUni.testMacs = HUni_test;
            databaseUni.trainingLabels = database.trainingLabels;
            databaseUni.testLabels = database.testLabels;
    
            [ results_kNNUni ] = knn_baseline( databaseUni , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'Uniform/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNUni');

            pos_error(iter,5) = mean(results_kNNUni.error(:,2));
            floor_rate(iter,5) = mean((results_kNNUni.error(:,4)==0).*(results_kNNUni.error(:,5)==0))*100;

            %% ELMInput initializer
            disp('----------------- ELMInput Initialization ---------------------');
            time02 = cputime;
            [BInp_train,ZHInp_train,HInp_train,RMSEInp_train, WInp_train]=ae_elm_encoder_w(trainingMacs,AF,neurons,Winp);   % Autoencoder
            timeCompTrainIMP = cputime-time02;
            time_compression_train(iter,1) = timeCompTrainIMP;

            time02 = cputime;
            [BInp_test,ZHInp_test,HInp_test,RMSEInp_test]=ae_elm_encoder_w(testMacs,AF,neurons,WInp_train);   % Autoencoder
            timeCompTestIMP = cputime-time02;
            time_compression_test(iter,1) = timeCompTestIMP;

            databaseInp.trainingMacs = HInp_train;
            databaseInp.testMacs = HInp_test;
            databaseInp.trainingLabels = database.trainingLabels;
            databaseInp.testLabels = database.testLabels;

            % Saving the compressed dataset
            save(strcat("compressed_databases/", databases{idx}, "_inp.mat"), "databaseInp", "-mat");
    
            [ results_kNNInp ] = knn_baseline( databaseInp , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'Input/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNInp');

            pos_error(iter,1) = mean(results_kNNInp.error(:,2));
            floor_rate(iter,1) = mean((results_kNNInp.error(:,4)==0).*(results_kNNInp.error(:,5)==0))*100;
            time_positioning(iter,1) = results_kNNInp.timefull(3);
            %% FID
            disp('----------------- FID Initialization ---------------------');
            [B_train,ZH_train,H_train,RMSE_train,W_train]=ae_elm_encoder_w(trainingMacs,AF,neurons,W);   % Autoencoder
            [B_test,ZH_test,H_test,RMSE_test,W_test]=ae_elm_encoder_w(testMacs,AF,neurons,W_train);   % Autoencoder
    
            databaseFID.trainingMacs = H_train;
            databaseFID.testMacs = H_test;
            databaseFID.trainingLabels = database.trainingLabels;
            databaseFID.testLabels = database.testLabels;
    
            [ results_kNNFID ] = knn_baseline( databaseFID , k(idx), representation{idx}, ['distance_' distance{idx}] );

            result_dir = ['Results/' databases{idx} '/' 'FID/' strcat(AF,'_',int2str(iter)) '/'];
            if ~exist(result_dir, 'dir')
                mkdir(result_dir)
            end
            save([result_dir strcat('KNN_', representation{idx},'_k_',int2str(k(idx)), '.mat') ], 'results_kNNFID');

            pos_error(iter,7) = mean(results_kNNFID.error(:,2));
            floor_rate(iter,7) = mean((results_kNNFID.error(:,4)==0).*(results_kNNFID.error(:,5)==0))*100;
        end
        % Orthogonal
        report(row_count,1) = mean(pos_error(:,1));
        report(row_count,2) = std(pos_error(:,1));
        report(row_count,3) = mean(floor_rate(:,1));
        report(row_count,4) = std(floor_rate(:,1));
        % He
        report(row_count,5) = mean(pos_error(:,2));
        report(row_count,6) = std(pos_error(:,2));
        report(row_count,7) = mean(floor_rate(:,2));
        report(row_count,8) = std(floor_rate(:,2));
        % Glorot
        report(row_count,9) = mean(pos_error(:,3));
        report(row_count,10) = std(pos_error(:,3));
        report(row_count,11) = mean(floor_rate(:,3));
        report(row_count,12) = std(floor_rate(:,3));
        % Gaussian
        report(row_count,13) = mean(pos_error(:,4));
        report(row_count,14) = std(pos_error(:,4));
        report(row_count,15) = mean(floor_rate(:,4));
        report(row_count,16) = std(floor_rate(:,4));
        % Uniform
        report(row_count,17) = mean(pos_error(:,5));
        report(row_count,18) = std(pos_error(:,5));
        report(row_count,19) = mean(floor_rate(:,5));
        report(row_count,20) = std(floor_rate(:,5));
         % Input
        report(row_count,1) = mean(pos_error(:,1));
        report(row_count,2) = std(pos_error(:,1));
        report(row_count,3) = mean(floor_rate(:,1));
        report(row_count,4) = std(floor_rate(:,1));
        report(row_count,5) = mean(time_initialization(:,1));
        report(row_count,6) = mean(time_compression_train(:,1));
        report(row_count,7) = mean(time_compression_test(:,1));
        report(row_count,8) = mean(time_positioning(:,1));
        % FID
        report(row_count,25) = mean(pos_error(:,7));
        report(row_count,26) = std(pos_error(:,7));
        report(row_count,27) = mean(floor_rate(:,7));
        report(row_count,28) = std(floor_rate(:,7));
        disp('');
    end
    disp('end');
end
save('ELM_FullDataset_Report.mat', "report");
