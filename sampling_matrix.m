clc;clear all
m=1000;
n=1000;
pp=5;
%A=randn(50,5)*randn(5,50);
for i=1:m
for j=1:n
A(i,j)=1/((i^pp+j^pp)^(1/pp));
end
end
% A=hilb(m);
rr=20;
% r=5;
[U,~,V]=svds(A,rr);
%top leverage scores
for i=1:m
    lr(i)=norm(U(i,:))^2;
end
for i=1:n
    lc(i)=norm(V(i,:))^2;
end

for r=1:rr
[~,r_1]=maxk(lr,r);
[~,r_2]=maxk(lc,r);
C=A(:,r_2);
R=A(r_1,:);
UU=pinv(C)*A*pinv(R);
E=A-C*UU*R;
E_1(r)=norm(E(:),'fro');
fprintf('Error for top leveraeg scores %d\n',norm(E(:)))
%leverage scores sampling
lrr=lr.^2;
lcc=lc.^2;
prob_1 = lrr / sum(lrr);
prob_2 = lcc / sum(lcc);
idx = randsample(m, r, true, prob_1);
idy = randsample(n, r, true, prob_2);
idx = unique(idx); % eliminate duplicates
idy = unique(idy); % eliminate duplicates
C = A(:, idy);
R = A(idx,:);
UU=pinv(C)*A*pinv(R);
E=A-C*UU*R;
E_2(r)=norm(E(:),'fro');
fprintf('Error of leverage score probability distribution %d\n',norm(E(:)))
%uniform sampling
r_1 = randsample(m,r);
r_2 = randsample(n,r);
C=A(:,r_2,:);
R=A(r_1,:,:);
UU=pinv(C)*A*pinv(R);
E=A-C*UU*R;
E_3(r)=norm(E(:),'fro');
%norm(E(:))
fprintf('Error of uniform samling %d\n',norm(E(:)))

[icol, irow, M]  = cur_deim(A, r);
C=A(:,icol,:);
R=A(irow,:,:);
UU=pinv(C)*A*pinv(R);
E=A-C*UU*R;
E_4(r)=norm(E(:),'fro');
fprintf('Error of the DEIM %d\n',norm(E(:)))
fprintf('Tubal rank %d\n',r)


[Fcol,column_ix] = bestcolumn_2(A,r);
[Fcol,row_ix] = bestcolumn_2(A',r);
C=A(:,column_ix,:);
R=A(row_ix,:,:);
UU=pinv(C)*A*pinv(R);
E=A-C*UU*R;
E_5(r)=norm(E(:),'fro');
fprintf('Error of the Best Sampling %d\n',norm(E(:)))
fprintf('Tubal rank %d\n',r)

end

rr=1:r;
semilogy(rr,E_1)
hold on
semilogy(rr,E_2)
hold on
semilogy(rr,E_3)
hold on
semilogy(rr,E_4)

semilogy(rr,E_5)
legend('Top leverage scores','Leverage scores sampling','Uniform sampling without replacement','DEIM','Best sampling');

function [icol, irow, M]  = cur_deim(A, k)

%CUR_DEIM  DEIM incurred CUR decomposition
% function [icol, irow, M] = cur_deim(A, k)
% icol contains the selected column indices 
% irow contains the selected row indices 
% M is the middle matrix of the CUR approximation 

% C = A(:,icol);  R = A(irow,:)
%
% Reference: Embree and Sorensen, 2016
% 
% (C) Perfect Gidisu, Michiel Hochstenbach 2020

[U, ~, V] = svds(A,k);

irow=zeros(1,k);
icol=zeros(1,k);
for j = 1:k
  [~, irow(j)] = max(abs(U(:,j)));
  [~, icol(j)] = max(abs(V(:,j)));
  if j<k
   U(:,j+1) = U(:,j+1) - U(:,1:j) * (U(irow(1:j),1:j) \ U(irow(1:j),j+1));
   V(:,j+1) = V(:,j+1) - V(:,1:j) * (V(icol(1:j),1:j) \ V(icol(1:j),j+1));
  end
end
C=A(:,icol);
R=A(irow,:);
M = pinv(C)*A*pinv(R);
end