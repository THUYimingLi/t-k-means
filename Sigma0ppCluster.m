function [ idx, mu, itr ] = Sigma0ppCluster( varargin )

if nargin == 2
    dataSet = varargin{1};
    k = varargin{2};
else
    return;
end


%Initialize variables
x = dataSet(:,:)';
p = size(x,1);
n = size(x,2);
g = k;
mu = zeros(p,g);

tau = zeros(g,n);
D = zeros(g,n);
c = (1+p)/2;
Q = 0;


%Initialize \mu
mu(:,1) = x(:,randperm(n,1));

for i = 2:g
    distance_matrix = zeros(i-1,n);
    for j = 1:n
        for m = 1:i-1
            distance_matrix(m,j) = sum((x(:,j)-mu(:,m)).^2);
        end
    end
    index = Roulettemethod(distance_matrix');
    mu(:,i) = x(:,index);
    clear distance_matrix;
end

mu = mu + 0.01;


count = 0;
for itr = 1:200
    %E-step
    for i = 1:g
        for j = 1:n
            D(i,j) = (x(:,j)-mu(:,i))'*(x(:,j)-mu(:,i));
        end
    end
    temp = 1./(D).^(c*ones(g,n));%g*n
    for i = 1:g
        for j = 1:n
            tau(i,j) = temp(i,j)/sum(temp(:,j));
        end
    end
    [~, index] = max(tau);
    r = zeros(g,n);
    for j = 1:n
        r(index(j),j) = 1;
    end
    
    
    %Optimization
    Q_old = Q;
    Q = 0;
    for i = 1:g
        for j = 1:n
            q = tau(i,j)*log(1/D(i,j)^c);
            if isnan(q) == 1
                Q = Q + 0;
            else
                Q = Q + q;
            end
        end
    end
    
    %M-step    
    for i = 1:g
        frac_m0 = zeros(p,1);
        frac_m1 = 0;
        for j = 1:n
            frac_m0 = frac_m0 + tau(i,j)*x(:,j);
            frac_m1 = frac_m1 + tau(i,j);
        end
        mu(:,i) = frac_m0/frac_m1;
    end
    mu(isnan(mu)) = 0; 

    %Print information
    disp(['itr£º',num2str(itr),'£¬Q£º',num2str(Q)]);
    if isnan(Q) == 1
        break;
    end
    if abs(Q_old-Q)<0.1
        count = count + 1;
    else
        count = 0;
    end
    if count > 3
        break;
    end
end


%Clustering
[~,index] = max(tau);
idx = index';
mu = mu';

end

function [index] = Roulettemethod(distance_matrix)

% Find shortest distance between one sample and its closest cluster centroid
[min_distance,~] = min(distance_matrix,[],2);

% Normalize for further operations
min_distance = min_distance ./ sum(min_distance);

% Construct roulette according to min_distance
temp_roulette = zeros(size(distance_matrix,1),1);
for i = 1:size(distance_matrix,1)
    temp_roulette(i,1) = sum(min_distance(1:i,:));
end

% Generate a random number for selection
temp_rand = rand();

% Find the corresponding index
for i = 1:size(temp_roulette,1)
    if((i == 1) && temp_roulette(i,1) > temp_rand)
        index = 1;
    elseif((temp_roulette(i,1) > temp_rand) && (temp_roulette(i-1,1) < temp_rand))
        index = i;
    end
end
end