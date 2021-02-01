function [ idx, mu, itr ] = kmedianCluster( varargin )

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
Q = 0;

tau = zeros(g,n);
D = zeros(g,n);


%Initialize parameters
median = randperm(n,g);


count = 0;
for itr = 1:200
    %E-step
    for i = 1:g
        D(i,:) = sum((x - x(:,median(1,i))*ones(1,n)).^2);
    end
    [~, index] = min(D);
    r = zeros(g,n);
    for j = 1:n
        r(index(j),j) = 1;
    end
    
    
    %Optimization
    Q_old = Q;
    Q = 0;
    for i = 1:g
        for j = 1:n
            Q = Q - r(i,j)*(x(:,j)-x(:,median(1,i)))'*(x(:,j)-x(:,median(1,i)));
        end
    end
    
	%M-step  
    for i = 1:g
        index = find(r(i,:)==1);
        rho = zeros(1,sum(r(i,:),2));
        for j = 1:sum(r(i,:),2)
            rho(1,j) = sum(sum(abs(x(:,index)-x(:,index(1,j)))));
        end
        [~,num] = min(rho');
        median(1,i) = index(1,num);
    end

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
[~,index] = min(D);
idx = index';
mu = x(:,median(1,:))';

end

