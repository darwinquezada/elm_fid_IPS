clear all; close all; clc;
load('databases/TUT1.mat');
addpath('positioning_functions','elm_functions');
conf = setConfig();
   
%% 1.- Transform data
database = get_data_representation(database, conf.dataRepresentation, ['distance_' conf.distFunction]);
%% 2.- Eliminar columnas con cero
count = 1;
for i=1:size(database.trainingMacs,2)
    if sum(database.trainingMacs(:,i))>0
        non_cero_ap(:,count) = database.trainingMacs(:,i);
        count = count + 1;
    end
end

dbtransformation = non_cero_ap;
dbtransformation(dbtransformation==0)=nan;

[maximum, index] = max(dbtransformation);

accepted_rss = 5; %N
accepted_cons_null_values = 3;
accepted_cut_points = 4;
cell_cut_points = {};
verf_cons_null = zeros(accepted_cons_null_values+1,1);
cluster = 1;

for i=1:size(dbtransformation,2)
    count_rss = 0;
    count_nan = 0;
    pos = 1;
    pos_cell = 1;
    cons_nan_condition = 0;
    for j=1:size(dbtransformation,1)
        if count_nan >= accepted_cons_null_values
            count_rss = 0;
            pos = 1;
            mcut_point = 0;
            if cons_nan_condition < 1
                pos_cell = pos_cell + 1;
                cluster = cluster + 1;
            end
            cons_nan_condition = cons_nan_condition + 1;
        end
        
        if ~isnan(dbtransformation(j,i))
            count_rss = count_rss + 1;
            count_nan = 0;
            for k = 1:size(dbtransformation,2)
                if ~isnan(dbtransformation(j,k))
                    if i ~= k
                        mcut_point(pos,1) = i; % AP de analisis
                        mcut_point(pos,2) = j; % Sample (row)
                        mcut_point(pos,3) = dbtransformation(j,i); % Sample RSS (row)
                        mcut_point(pos,4) = k; % AP de corte
                        mcut_point(pos,5) = dbtransformation(j,k); % AP de corte
                        pos = pos + 1;
                    end
                end
            end
        else
            count_nan = count_nan + 1;
            for k = 1:size(dbtransformation,2)
                if ~isnan(dbtransformation(i,k))
                    if i ~= k
                        mcut_point(pos,1) = i; % AP de analisis
                        mcut_point(pos,2) = j; % Sample (row)
                        mcut_point(pos,3) = dbtransformation(j,i); % Sample RSS (row)
                        mcut_point(pos,4) = k; % AP de corte
                        mcut_point(pos,5) = dbtransformation(j,k); % AP de corte
                        pos = pos + 1;
                    end
                end
            end
        end
        
        if count_rss >= accepted_rss && count_nan <= accepted_cons_null_values
            cons_nan_condition = 0;
            if j<(size(dbtransformation,1)-accepted_cons_null_values) && j>1
                temp = j-1;
                for k = 1:accepted_cons_null_values+1
                    if isnan(dbtransformation(temp+k,i))
                        verf_cons_null(k) = 1;
                    end
                end
                if (size(find(verf_cons_null==1),2)>accepted_cons_null_values)
                    mcut_point = 0;
                end
            end

            for l =1:size(mcut_point,1)
                mcut_point(l,6) = cluster;
            end
            if mcut_point ~= 0
                cell_cut_points{pos_cell,i} = mcut_point;
            end
            verf_cons_null = 0;
        end
    end
end

cluster = 1;
for i = 1:size(cell_cut_points,2)
    for j = 1:size(cell_cut_points,1)
        if size(cell_cut_points{j,i},1)~=0
            unique_samples = unique(cell_cut_points{j,i}(:,2));
            % Contar ocurrencias por AP
            srt = cell_cut_points{j,i}(:,4);
            unique_apc = unique(cell_cut_points{j,i}(:,4));
            GC = zeros(size(unique(srt),1),1); % Group count
            GR = zeros(size(unique(srt),1),1); % Group row
            GC_APC = zeros(size(unique(srt),1),1); % Group count

            for k = 1:size(GC,1)
                GR(k) = unique_apc(k);
                %before = 0;
                count = 1;
                init = 1;
                max_rss = 0;
                for x = 1:size(srt,1)
                    if GR(k) == srt(x)
                        GC(k) = count;
                        count = count+1;
                        
                        if init == 1
                            max_rss = cell_cut_points{j,i}(x,5);
                            init = 0;
                        end
                        if max_rss > cell_cut_points{j,i}(x,5)
                            GC_APC(k) = max_rss;
                            %max_rss = max_rss;
                        else
                            GC_APC(k) = cell_cut_points{j,i}(x,5);
                            max_rss = cell_cut_points{j,i}(x,5);
                        end
                        %before = cell_cut_points{j,i}(x,5);
                    end
                end
            end

            % Maximum by APC

%             for k = 1:size(GC,1)
%                 before = 0;
%                 for x = 1:size(srt,1)
%                     if GR(k) == srt(x)
%                         if k == 1
%                             maximum = cell_cut_points{j,i}(x,5);
%                         end
%                         if maximum > before
%                             GC_APC(k) = cell_cut_points{j,i}(x,5);
%                             maximum = cell_cut_points{j,i}(x,5);
%                         else
%                             GC_APC(k) = before;
%                             maximum = before;
%                         end
%                         before = cell_cut_points{j,i}(x,5);
%                     end
%                 end
%             end

            aver_apc = mean(GC);
            aver_apc_rss = GC_APC;
            max_apc_rss = max(GC_APC);
        end
    end
end 

%% 3.- Determinar muestras dispersas

%% 4.- Promedio entre las columnas con mayor grado de afinidad

%% 5.- Reducir el número de columans al número de neuronas

%% 6.- Multiplicar las columas por los coeficientes del PCA (Buscar mejores métodos)
