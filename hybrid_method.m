clc;clear all;close all
A=double(imread('peppers.bmp'));
[m,n,p]=size(A);
rr=50;
r=rr;
ll=5;

[U,~,V]=tensor_t_svd(A,rr-ll);
% 
% irow=zeros(1,r);
% icol=zeros(1,r);
% for j = 1:r
%     for i=1:m
%     n_1(i)=norm(squeeze(U(i,j,:)),'fro');
%     end
%     for i=1:n
%     n_2(i)=norm(squeeze(V(i,j,:)),'fro');
%     end
%   [~, irow(j)] = max(n_1);
%   [~, icol(j)] = max(n_2);
%   if j<r
%    U(:,j+1,:) = U(:,j+1,:) - tprod(U(:,1:j,:),tprod(t_pinv(U(irow(1:j),1:j,:)),U(irow(1:j),j+1,:))); 
%    V(:,j+1,:) = V(:,j+1,:) - tprod(V(:,1:j,:),tprod(t_pinv(V(icol(1:j),1:j,:)),V(icol(1:j),j+1,:)));
%   end
% end

[irow,icol]=tdeim(A,r);

for i=1:m
l_1(i)=norm(squeeze(U(i,:,:)),'fro');
end
% l_1=sort(l_1,'descend');
% [nn_1,nn_2]=maxk(l_1,10)
l_1=sort(l_1);
l_1(irow)=0;
[~,pp_1]=maxk(l_1,ll);
l_1=union(pp_1,irow);

for i=1:n
l_2(i)=norm(squeeze(V(i,:,:)),'fro');
end
% l_2=sort(l_2);
% [mm_1,mm_2]=maxk(l_1,10);
l_2=sort(l_2);
l_2(icol)=0;
[~,pp_2]=maxk(l_2,ll);
l_2=union(pp_2,icol);

C=A(:,l_2,:);
R=A(l_1,:,:);
M = tprod(tprod(t_pinv(C),A),t_pinv(R));
E=A-tprod(tprod(C,M),R);
%norm(E(:))
fprintf('Error of the hybrid DEIM %d\n',norm(E(:)))
%%
 [U,~,V]=tensor_t_svd(A,rr);

 [irow,icol]=tdeim(A,r);

C=A(:,icol,:);
R=A(irow,:,:);
M = tprod(tprod(t_pinv(C),A),t_pinv(R));
E=A-tprod(tprod(C,M),R);
%norm(E(:))
fprintf('Error of the DEIM %d\n',norm(E(:)))