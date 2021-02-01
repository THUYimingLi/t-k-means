function [ idx, mu, itr ] = gmmCluster( varargin )
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
Pi(1:g,1) = 1/g;
for i = 1:g
    Sigma(:,:,i) = eye(p);
end

tau = zeros(g,n);
D = zeros(g,n);
Q = 0;
dQ = 0;


%Initialize parameters
low = min(x,[],2)*ones(1,g);
up = max(x,[],2)*ones(1,g);
mu(:,:) = unifrnd(low,up);


count = 0;
for itr = 1:200
    %E-step
    temp = zeros(g,n);
    for i = 1:g
        for j = 1:n
            D(i,j) = (x(:,j)-mu(:,i))'/Sigma(:,:,i)*(x(:,j)-mu(:,i));
            temp(i,j) = Pi(i,1)/det(Sigma(:,:,i))^0.5/exp(-0.5*D(i,j));
        end
    end
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
            q = tau(i,j)*log(Pi(i,1)/det(Sigma(:,:,i))^0.5/exp(-0.5*D(i,j)));
            if isnan(q) == 1
                Q = Q + 0;
            else
                Q = Q + q;
            end
        end
    end
    
    %M-step
    frac_a0 = 0;
    frac_a1 = 0;
    for i = 1:g
        Pi(i,1) = sum(tau(i,:))/n;
        frac_m0 = zeros(p,1);
        frac_m1 = 0;
        for j = 1:n
            frac_m0 = frac_m0 + tau(i,j)*x(:,j);
            frac_m1 = frac_m1 + tau(i,j);
        end
        mu(:,i) = frac_m0/frac_m1;
        for j = 1:n
            frac_a0 = frac_a0 + tau(i,j)*(x(:,j)-mu(:,i))*(x(:,j)-mu(:,i))';
            frac_a1 = frac_a1 + tau(i,j);
        end
        if det(frac_a0/frac_a1)>0.0000000001
            Sigma(:,:,i) = frac_a0/frac_a1;
        else
            Sigma(:,:,i) = frac_a0/frac_a1 + 0.0000000001*eye(p);
        end
    end
    mu(isnan(mu)) = 0;

    %Print information
    disp(['itr£º',num2str(itr),'£¬Q£º',num2str(Q)]);
    if isnan(Q) == 1
        break;
    end
    if abs(Q_old-Q)<0.1 || abs(dQ - abs(Q_old-Q))<0.1
        count = count + 1;
    else
        count = 0;
    end
    dQ=abs(Q_old-Q);
    if count > 3
        break;
    end
end

%Clustering
[~,index] = max(tau);
idx = index';
mu = mu';

end

