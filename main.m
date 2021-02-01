warning('off');
addpath('.\libs');


%Specify experiment settings
data_version = 'dim064';
method_count = 2; %Evaluate two methods simultaneously
exp_count = 10; %Repeat the experiment for 10 times

%Obtain the name of datasets
dataNames = who('-file',['datamat/',data_version ,'.mat']);


for exp_itr = 1:exp_count
    clearvars -EXCEPT exp_itr data_version dataNames method_count exp_count;
    
    %Load all datasets
    load(['datamat/', data_version, '.mat']);
 
    for n = 1:size(dataNames,1)
        %obtain dataset
        dataSet = {};
        dataName = dataNames{n};
        eval(['dataSet = ' dataName ';']);

        %Obtain #clusters k
        centroids = dataSet{2};
        k = size(centroids,1);

        %Obtain samples
        [samples,~] = mapminmax((zscore(dataSet{1}(:,1:end-1)))');
        samples = samples';
        labels = dataSet{1}(:,end);
        noise = unifrnd(-1,1,ceil(size(samples,1)*0.2),size(samples,2));
        trainset = [samples;noise];

        %Clustering
        idx = cell(method_count,1);
        mu = cell(method_count,1);
        itr = zeros(method_count,1);
        time = zeros(method_count,1);
        
		%---------------------------One may modify this part to evaluate other methods---------------------------
        tic
        [idx{1,1},mu{1,1},itr(1,1)] = kmeansCluster(trainset,k);
        idx{1,1} = idx{1,1}(1:end-size(noise,1),:);
        time(1,1) = toc;
        tic
        [idx{2,1},mu{2,1},itr(2,1)] = Sigma0Cluster(trainset,k);
        idx{2,1} = idx{2,1}(1:end-size(noise,1),:);
        time(2,1) = toc;
		%---------------------------One may modify this part to evaluate other methods---------------------------
		
        %Visualization
        figure('NumberTitle', 'off', 'Name', dataName);
        for i = 1:method_count
            subplot(1,2,i);
            scatter(samples(:,1),samples(:,2),10,'x');
            hold on;
            voronoi(mu{i,1}(:,1),mu{i,1}(:,2),'k');
            hold on;
            scatter(mu{i,1}(:,1),mu{i,1}(:,2),'y*');
            hold on;
            axis square;
            title(['ARI:', num2str(RandIndex(labels,idx{i,1}))]);
        end

        %Calculate Evaluation Metric
            %ARI
            ARI = zeros(method_count,1);
            for z = 1:method_count
                ARI(z,1) = RandIndex(labels,idx{z,1});
%                 ARI(z,1) = 0;
            end
            %NMI
            nmi = zeros(method_count,1);
            for z = 1:method_count
                nmi(z,1) = NMI(labels,idx{z,1});
%                 NMI(z,1) = 0;
            end
            %WSS BSS
            WSS = zeros(method_count,1);
            BSS = zeros(method_count,1);
            for z = 1:method_count
                for x = 1:size(samples,1)
                    for y = 1:size(samples,2)
                        WSS(z,1) = WSS(z,1)+(samples(x,y)-mu{z,1}(idx{z,1}(x,1),y))^2;
                    end
                end
                WSS(z,1) = WSS(z,1)/size(idx{z,1}(:,1),1);
            end
            Q = mean(samples);
            for z = 1:method_count
                for x = 1:k
                    for y = 1:size(samples,2)
                        BSS(z,1) = BSS(z,1)+size(find(idx{z,1}(:,1)==x),1)*(Q(1,y)-mu{z,1}(x,y))^2;
                    end
                end
                BSS(z,1) = BSS(z,1)/size(idx{z,1}(:,1),1);
            end
            WB = WSS(:,:)./BSS(:,:);
            %DbDunn
            cintra = zeros(k,1,method_count);
            cinter = zeros(k,k,method_count);
            DB = zeros(method_count,1);
            D = zeros(method_count,1);
            for z = 1:method_count
                for x = 1:k
                    cintra(x,1,z) = mean(sqrt(sum(abs(samples(find(idx{z,1}==x),:)-ones(size(find(idx{z,1}==x),1),1)*mean(samples(find(idx{z,1}==x),:))).^2,2)));
                    for y = 1:k
                        cinter(x,y,z) = sqrt(sum(abs(mean(samples(find(idx{z,1}==x),:))-mean(samples(find(idx{z,1}==y),:))).^2,2));
                    end
                end
            end
            for z = 1:method_count
                [DB(z,1),D(z,1)] = DbDunn(cintra(:,:,z),cinter(:,:,z),k);
                if isnan(D(z,1))
                    D(z,1) = 3;
                end
            end
            %MSE

                MSE = zeros(method_count,1);
                for z = 1:method_count
                    for x = 1:k
                        MSE(z,1) = MSE(z,1)+sum(sum(abs(samples(find(idx{z,1}==x),:)-ones(size(find(idx{z,1}==x),1),1)*mean(samples(find(idx{z,1}==x),:))).^2,2));
                    end
                    MSE(z,1) = MSE(z,1)/size(idx{z,1}(:,1),1);
                end

            save(['results/', dataName,'_', num2str(exp_itr),'.mat'],'MSE','WB','DB','D','ARI','nmi','itr','time');

    end
end

clearvars -EXCEPT data_version dataNames method_count exp_count;
index_count = 8;
table = cell(method_count,size(dataNames,1),index_count);
for exp_itr = 1:exp_count
    for data_itr = 1:size(dataNames,1)
        load(['results/',dataNames{data_itr},'_',num2str(exp_itr),'.mat']);
        for method_itr = 1:method_count
            table{method_itr,data_itr,1} = [table{method_itr,data_itr,1} MSE(method_itr,1)];
            table{method_itr,data_itr,2} = [table{method_itr,data_itr,2} WB(method_itr,1)];
            table{method_itr,data_itr,3} = [table{method_itr,data_itr,3} DB(method_itr,1)];
            if isnan(D(method_itr,1))
                D(method_itr,1) = 0;
            end
            table{method_itr,data_itr,4} = [table{method_itr,data_itr,4} D(method_itr,1)];
            table{method_itr,data_itr,5} = [table{method_itr,data_itr,5} ARI(method_itr,1)];
            table{method_itr,data_itr,6} = [table{method_itr,data_itr,6} nmi(method_itr,1)];
            table{method_itr,data_itr,7} = [table{method_itr,data_itr,7} itr(method_itr,1)];
            table{method_itr,data_itr,8} = [table{method_itr,data_itr,8} time(method_itr,1)];
        end
    end
end
for data_itr = 1:size(dataNames,1)
    for method_itr = 1:method_count
        for index = 1:3
            table{method_itr,data_itr,index} = [min(table{method_itr,data_itr,index}) mean(table{method_itr,data_itr,index}) std(table{method_itr,data_itr,index})];
        end
        table{method_itr,data_itr,4} = [max(table{method_itr,data_itr,4}) mean(table{method_itr,data_itr,4}) std(table{method_itr,data_itr,4})];
        table{method_itr,data_itr,5} = [max(table{method_itr,data_itr,5}) mean(table{method_itr,data_itr,5}) std(table{method_itr,data_itr,5})];
        table{method_itr,data_itr,6} = [max(table{method_itr,data_itr,6}) mean(table{method_itr,data_itr,6}) std(table{method_itr,data_itr,6})];
        table{method_itr,data_itr,7} = [min(table{method_itr,data_itr,7}) mean(table{method_itr,data_itr,7}) std(table{method_itr,data_itr,7})];
        table{method_itr,data_itr,8} = [min(table{method_itr,data_itr,8}) mean(table{method_itr,data_itr,8}) std(table{method_itr,data_itr,8})];
    end
end
save('results/table.mat','table');
figure;
excel = cell(index_count,3);
for i = 1:index_count
    for j = 1:3
        excel{i,j}=zeros(method_count,size(dataNames,1));
    end
end
for index = 1:index_count
    for data_itr = 1:size(dataNames,1)
        subplot(index_count,size(dataNames,1),size(dataNames,1)*(index-1)+ data_itr);
        for method_itr = 1:method_count
            best_value = table{method_itr,data_itr,index}(1,1);
            mean_value = table{method_itr,data_itr,index}(1,2);
            std_value = table{method_itr,data_itr,index}(1,3);
            errorbar(method_itr,mean_value,std_value,'o');
            hold on;
            scatter(method_itr,best_value,'*');
            hold on;
            excel{index,1}(method_itr,data_itr) = best_value;
            excel{index,2}(method_itr,data_itr) = mean_value;
            excel{index,3}(method_itr,data_itr) = std_value;
        end
%         set(gca,'xtick',[1 2 3 4 5 6 7])
%         set(gca,'xticklabel',{'k-means','k-means++','k-Medoid','GMM','TMM','SigmaAlpha','Sigma0'})
    end
end
save('results/excel.mat','excel');