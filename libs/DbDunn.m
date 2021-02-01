function [DB, Dunn] = DbDunn(cintra, cinter, k)
% Davies-Bouldin index
  R = zeros(k);
  dbs=zeros(1,k);
  for i = 1:k
    for j = i+1:k
      if cinter(i,j) == 0 
         R(i,j) = 0;
      else
         R(i,j) = (cintra(i) + cintra(j))/cinter(i,j);
      end
    end
    dbs(i) = max(R(i,:));
  end
  DB = mean(dbs(1:k-1));
  
  % Dunn index
  dbs = max(cintra);
  R = cinter/dbs;
  for i = 1:k-1
     S = R(i,i+1:k);
     dbs(i) = min(S);
  end
  Dunn = min(dbs); 
 