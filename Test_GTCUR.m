clc;clear all;close all
%A=randn(20,20,20);
%B=randn(30,20,20);
A=double(imread('lena.bmp'));
B=double(imread('peppers.bmp'));

r=20;
[U,V,X,C,S]=gsvd_tube(A,B);
U_1=U(:,1:r,:);
V_1=V(:,1:r,:);
X_1=X(:,1:r,:);
C_1=C(1:r,1:r,:);
S_1=S(1:r,1:r,:);
A_1=t_prod(t_prod(U_1,C_1),t_trans(X_1));
B_1=t_prod(t_prod(V_1,S_1),t_trans(X_1));
EE_1=A-A_1;
norm(EE_1(:))
EE_2=B-B_1;
norm(EE_2(:))


subplot(2,2,1)
imshow(uint8(A_1))
subplot(2,2,2)
imshow(uint8(B_1))


[UU,SS,VV]=tensor_t_svd(A,r);
subplot(2,2,3)
imshow(uint8(t_prod(t_prod(UU,SS),t_trans(VV))))
[UU_1,SS_1,VV_1]=tensor_t_svd(B,r);
subplot(2,2,4)
imshow(uint8(t_prod(t_prod(UU_1,SS_1),t_trans(VV_1))))
A_1=(t_prod(t_prod(UU,SS),t_trans(VV)));
B_1=(t_prod(t_prod(UU_1,SS_1),t_trans(VV_1)));
EE_1=A-A_1;
norm(EE_1(:))
EE_2=B-B_1;
norm(EE_2(:))

%%
% r=60;
[irow_1,~]=tdeim(U_1,r);
[irow_2,~]=tdeim(V_1,r);
[~,icol]=tdeim(X,r);
R_1=A(irow_1,:,:);
R_2=B(irow_2,:,:);
C_1=A(:,icol,:);
C_2=B(:,icol,:);
M_1=t_prod(t_prod(t_pinv(C_1),A),t_pinv(R_1));
M_2=t_prod(t_prod(t_pinv(C_2),B),t_pinv(R_2));
A_ap=t_prod(t_prod(C_1,M_1),R_1);
ER=A-A_ap;
norm(ER(:))

B_ap=t_prod(t_prod(C_2,M_2),R_1);
ER=B-B_ap;
norm(ER(:))

subplot(1,2,1)
imshow(uint8(A_ap))
subplot(1,2,2)
imshow(uint8(B_ap))

%%
[m,n]=size(A);
for i=1:m
    lr_1(i)=norm(squeeze(U_1(i,:,:)),'fro');
end

for i=1:m
    lr_2(i)=norm(squeeze(UU(i,:,:)),'fro');
end

figure(1)
%stem(lr,'*')
plot(lr_1,'o')
hold on
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
legend('Tubal leverage scores (horizontal)','TDEIM horizontal slice indices','Position',[0.44 0.92 0.15 0.0869]);
legend('Orientation','horizontal')
xlim([0 256])

figure(2)
%stem(lr,'*')
plot(lr_2,'o')
hold on
A=zeros(m,0);
A(irow_1,1)=1;
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
legend('Tubal leverage scores (horizontal)','TDEIM horizontal slice indices','Position',[0.44 0.92 0.15 0.0869]);
legend('Orientation','horizontal')
xlim([0 256])
%%

for i=1:n
    lc(i)=norm(squeeze(V_1(i,:,:)),'fro');
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