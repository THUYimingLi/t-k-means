function [ idx, mu, itr ] = kmeansCluster( varargin )

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
Pi = 1/g;
mu = zeros(p,g);

r = zeros(g,n);
tau = zeros(g,n);
Q = 0;
Q_old = 0;


%Initialize parameters
low = min(x,[],2)*ones(1,g);
up = max(x,[],2)*ones(1,g);
mu(:,:) = unifrnd(low,up);


count = 0;
for itr = 1:200
    %E-step
    for i = 1:g
        for j = 1:n
            tau(i,j) = (x(:,j)-mu(:,i))'*(x(:,j)-mu(:,i));
        end
    end
    [~, index] = min(tau);
    r = zeros(g,n);
    for j = 1:n
        r(index(j),j) = 1;
    end
    
    
    %Optimization
    Q_old = Q;
    Q = 0;
    for i = 1:g
        for j = 1:n
            Q = Q - r(i,j)*(x(:,j)-mu(:,i))'*(x(:,j)-mu(:,i));
        end
    end
    
    %M-step
    for i = 1:g
        mu(:,i) = mean(x(:,find(r(i,:)==1)),2);
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
[~,index] = min(tau);
idx = index';
mu = mu';

end

