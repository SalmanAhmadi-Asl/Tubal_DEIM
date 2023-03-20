clc;clear all;close all
m=300;
n=400;
p=300;
rr=15;
r=rr;
%A=tprod(randn(m,rr,p),randn(rr,n,p))+.1*randn(m,n,p);
% for i=1:p
% A(:,:,i)=hilb(m);
% end
pp=3;
for i=1:m
    for j=1:n
        for k=1:p
A(i,j,k)=1/((i^pp+j^pp+k^pp)^(1/pp));
        end
    end
end


[U,~,V]=tensor_t_svd(A,r);
%% tubal leverage scores (top ones)
for i=1:m
    lr(i)=norm(squeeze(U(i,:,:)),'fro');
end
for i=1:n
    lc(i)=norm(squeeze(V(i,:,:)),'fro');
end
%subplot(1,2,1)
figure(1)
%stem(lr,'*')
plot(lr,'o')
hold on
%subplot(1,2,2)
figure(2)
%stem(lc,'o')
plot(lc,'o')
hold on
[~,r_1]=maxk(lr,r);
[~,r_2]=maxk(lc,r);
C=A(:,r_2,:);
R=A(r_1,:,:);
UU=tprod(tprod(t_pinv(C),A),t_pinv(R));
E=A-tprod(tprod(C,UU),R);
fprintf('Error for top leveraeg scores %d\n',norm(E(:)))
%% sampling with  tubal leverage scores probability distribution
lrr=lr.^2;
lcc=lc.^2;
prob_1 = lrr / sum(lrr);
prob_2 = lcc / sum(lcc);
idx = randsample(m, r, true, prob_1);
idy = randsample(n, r, true, prob_2);
idx = unique(idx); % eliminate duplicates
idy = unique(idy); % eliminate duplicates
C = A(:, idy,:);
R = A(idx,:,:);
UU=tprod(tprod(t_pinv(C),A),t_pinv(R));
E=A-tprod(tprod(C,UU),R);
fprintf('Error of leverage score probability distribution %d\n',norm(E(:)))
%% uniform sampling
r_1 = randsample(m,r);
r_2 = randsample(n,r);
C=A(:,r_2,:);
R=A(r_1,:,:);
UU=tprod(tprod(t_pinv(C),A),t_pinv(R));
E=A-tprod(tprod(C,UU),R);
%norm(E(:))
fprintf('Error of uniform samling %d\n',norm(E(:)))
%% tubal DEIM
% clc;clear all
% m=50;
% n=50;
% p=50;
% r=5;
% A=tprod(randn(m,r,p),randn(r,n,p));
% [U,~,V]=tensor_t_svd(A,10);
irow=zeros(1,r);
icol=zeros(1,r);
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
C=A(:,icol,:);
R=A(irow,:,:);
M = tprod(tprod(t_pinv(C),A),t_pinv(R));
E=A-tprod(tprod(C,M),R);
%norm(E(:))
fprintf('Error of the DEIM %d\n',norm(E(:)))
%%
A=zeros(m,0);
A(irow,1)=1;
for i=1:m
    if A(i,1)==0
        A(i,1)=nan;
    end
end
%subplot(1,2,1)
figure(1)
H=stem(A)
set(H, 'Marker', 'none','LineWidth',.3)
xlabel('Horizontal slice (feature) index')
legend('TDEIM horizontal slice indices','Tubal leverage scores','Position',[0.44 0.92 0.15 0.0869]);
legend('Orientation','horizontal')
A=zeros(m,0);
A(icol,1)=1;
for i=1:m
    if A(i,1)==0
        A(i,1)=nan;
    end
end
%subplot(1,2,2)
figure(2)
H=stem(A)
set(H, 'Marker', 'none','LineWidth',.3)
xlabel('Lateral slice (feature) index')
legend('TDEIM lateral slice indices','Tubal leverage scores','Position',[0.44 0.92 0.15 0.0869]);
legend('Orientation','horizontal')
%irow-icol
% [icol, irow, M]  = cur_deim(A, r);
% warning off
% 
% function [icol, irow, M]  = cur_deim(A, r)
% 
% [U, ~, V] = tensor_t_svd(A,r);
% 
% irow=zeros(1,r);
% icol=zeros(1,r);
% for j = 1:1
%   [~, irow(j)] = max(norm(squeeze(U(:,j,:)),'fro'));
%   [~, icol(j)] = max(norm(squeeze(V(:,j,:)),'fro'));
%   if j<r
% %       tinv(U(irow(1:j),1:j,:))
%    U(:,j+1,:) = U(:,j+1,:) - tprod(U(:,1:j,:),tprod(t_pinv(U(irow(1:j),1:j,:)),U(irow(1:j),j+1,:))); 
%    V(:,j+1,:) = V(:,j+1,:) - tprod(V(:,1:j,:),tprod(t_pinv(V(icol(1:j),1:j,:)),V(icol(1:j),j+1,:)));
%    
%    
% %    V(:,1:j,:) * (V(icol(1:j),1:j,:) \ V(icol(1:j),j+1,:));
%   end
% end
% C=A(:,icol,:);
% R=A(irow,:,:);
% M = tprod(tprod(t_pinv(C),A),t_pinv(R));
% end