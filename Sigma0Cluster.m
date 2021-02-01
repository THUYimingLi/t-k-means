function [ idx, mu, itr ] = Sigma0Cluster( varargin )

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


%Initialize parameters
low = min(x,[],2)*ones(1,g); 
up = max(x,[],2)*ones(1,g);
mu(:,:) = unifrnd(low,up);



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

