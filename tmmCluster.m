function [ idx, mu, itr ] = tmmCluster( varargin )

if nargin == 2 || nargin == 3
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
nu(1:g,1) = 15;
Pi(1:g,1) = 1/g;
for i = 1:g
    Sigma(:,:,i) = eye(p);
end
if nargin == 3
    nu = varargin{3};
end

tau = zeros(g,n);
rho = zeros(g,n);
lrho = zeros(g,n);
D = zeros(g,n);
Q = 0;
dQ = 0;


%Initialize parameters
low = min(x,[],2)*ones(1,g);
up = max(x,[],2)*ones(1,g);
mu(:,:) = unifrnd(low,up);


count = 0;
for itr = 1:500
    %E-step
    temp = zeros(g,n);
    for i = 1:g
        for j = 1:n
            D(i,j) = (x(:,j)-mu(:,i))'/Sigma(:,:,i)*(x(:,j)-mu(:,i));
            temp(i,j) = Pi(i,1)*gamma((nu(i,1)+p)/2)/gamma(nu(i,1)/2)/det(Sigma(:,:,i))^0.5/(pi*nu(i,1))^(p/2)/(1+D(i,j)/nu(i,1))^((nu(i,1)+p)/2);
        end
    end
    for i = 1:g
        for j = 1:n
            tau(i,j) = temp(i,j)/sum(temp(:,j));
            rho(i,j) = (nu(i,1)+p)/(nu(i,1)+D(i,j));
            lrho(i,j) = log(rho(i,j)) - rho(i,j);
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
            q = tau(i,j)*log(Pi(i,1)*gamma((nu(i,1)+p)/2)/gamma(nu(i,1)/2)/det(Sigma(:,:,i))^0.5/(pi*nu(i,1))^(p/2)/(1+D(i,j)/nu(i,1))^((nu(i,1)+p)/2));
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
            frac_m0 = frac_m0 + tau(i,j)*x(:,j)/rho(i,j);
            frac_m1 = frac_m1 + tau(i,j)/rho(i,j);
        end
        mu(:,i) = frac_m0/frac_m1;
        for j = 1:n
            frac_a0 = frac_a0 + tau(i,j)*rho(i,j)*(x(:,j)-mu(:,i))*(x(:,j)-mu(:,i))';
            frac_a1 = frac_a1 + tau(i,j);
        end
        Sigma(:,:,i) = frac_a0/frac_a1 + 0.0001*eye(p);
    end
    mu(isnan(mu)) = 0;
    
    for i = 1:g
        temp = 0;
        for j = 1:n
            temp = temp + tau(i,j)*(log((nu(i,1)+p)/(nu(i,1)+(x(:,j)-mu(:,i))'/Sigma(:,:,i)*(x(:,j)-mu(:,i)))) - (nu(i,1)+p)/(nu(i,1)+(x(:,j)-mu(:,i))'/Sigma(:,:,i)*(x(:,j)-mu(:,i))));
        end
        star = temp/sum(tau(i,:),2);
        nu_old = nu(i,1);
        nu(i,1) = 1/(-1-star+log((nu(i,1)+p)/2)-psi((nu(i,1)+p)/2));

        delta = 10;
        while(delta>0)
            nu(i,1) = nu(i,1)+0.01;
            delta = -psi(nu(i,1))+log(nu(i,1))+1+star-log((nu_old+p)/2)+psi((nu_old+p)/2);
        end
    end
    
    
    if nargin == 3
        nu = varargin{3};
    end

    %Print information
    disp(['itr£º',num2str(itr),'£¬Q£º',num2str(Q),'£¬nu£º',num2str(sum(nu)/g)]);
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

