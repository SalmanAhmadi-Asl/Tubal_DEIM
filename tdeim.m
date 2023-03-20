function [irow,icol]=tdeim(A,r)
[U,~,V]=tensor_t_svd(A,r);
irow=zeros(1,r);
icol=zeros(1,r);
[m,n,~]=size(A);
for j = 1:r
    for i=1:m
    n_1(i)=norm(squeeze(U(i,j,:)),'fro');
    end
    for i=1:n
    n_2(i)=norm(squeeze(V(i,j,:)),'fro');
    end
  [~, irow(j)] = max(n_1);
  [~, icol(j)] = max(n_2);
  if j<r
   U(:,j+1,:) = U(:,j+1,:) - tprod(U(:,1:j,:),tprod(t_pinv(U(irow(1:j),1:j,:)),U(irow(1:j),j+1,:))); 
   V(:,j+1,:) = V(:,j+1,:) - tprod(V(:,1:j,:),tprod(t_pinv(V(icol(1:j),1:j,:)),V(icol(1:j),j+1,:)));
  end
end
