function [ idx, mu, itr ] = SigmaAlphaCluster( varargin )

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
nu = 15;
alpha = 0.01;
if nargin == 3
    nu = varargin{3};
end

tau = zeros(g,n);
rho = zeros(g,n);
lrho = zeros(g,n);
D = zeros(g,n);
c = (p+nu)/2;
Q = 0;
dQ = 0;


%Initialize parameters
low = min(x,[],2)*ones(1,g);
up = max(x,[],2)*ones(1,g);
mu(:,:) = unifrnd(low,up);



count = 0;
for itr = 1:500
    %E-step
    for i = 1:g
        for j = 1:n
            D(i,j) = (x(:,j)-mu(:,i))'*(x(:,j)-mu(:,i));
        end
    end
    temp = 1./(1+D/(nu*alpha)).^(c*ones(1,n));%g*n
    for i = 1:g
        for j = 1:n
            tau(i,j) = temp(i,j)/sum(temp(:,j));
            rho(i,j) = (nu+p)/(nu+D(i,j)/alpha);
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
            q = tau(i,j)*log(1/(1+D(i,j)/(nu*alpha))^c);
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
        frac_m0 = zeros(p,1);
        frac_m1 = 0;
        for j = 1:n
            frac_m0 = frac_m0 + tau(i,j)*x(:,j)/rho(i,j);
            frac_m1 = frac_m1 + tau(i,j)/rho(i,j);
        end
        mu(:,i) = frac_m0/frac_m1;
        for j = 1:n
            frac_a0 = frac_a0 + tau(i,j)*rho(i,j)*D(i,j);
            frac_a1 = frac_a1 + tau(i,j);
        end
    end
    alpha = frac_a0/frac_a1/p;
    mu(isnan(mu)) = 0;
    
    star = 0;
    for i = 1:g
        temp = 0;
        for j = 1:n
            temp = temp + tau(i,j)*(log((nu+p)/(nu+(x(:,j)-mu(:,i))'*(x(:,j)-mu(:,i))/alpha)) - (nu+p)/(nu+(x(:,j)-mu(:,i))'*(x(:,j)-mu(:,i))/alpha));
        end
        star = star + temp/sum(tau(i,:),2);
    end
    star = star/g;
    nu_old = nu;
    nu = 1/(-1-star+log((nu+p)/2)-psi((nu+p)/2));

    delta = 10;
    while(delta>0)
        nu = nu+0.01;
        delta = -psi(nu)+log(nu)+1+star-log((nu_old+p)/2)+psi((nu_old+p)/2);
    end
    
    if nargin == 3
        nu = varargin{3};
    end

    %Print information
    disp(['itr£º',num2str(itr),'£¬Q£º',num2str(Q),'£¬nu£º',num2str(nu),'£¬alpha£º',num2str(alpha)]);
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

