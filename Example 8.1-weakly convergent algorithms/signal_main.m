
function signal_main
clc
tic
% n is the original signal length
m =480;%1440; %% 240  480 720 960 1200
% k is number of observations to make
k =2048;%6144; %  1024 2048 3072 4096 5120
% number of spikes to put down
p=60; %180; % 30  60 90 120 150
f = zeros(k,1);
q = randperm(k);
f(q(1:p)) = sign(randn(p,1)); %%  original signal
% measurement matrix
R = randn(m,k);  
% orthonormalize rows
R = orth(R')';
b=R*f; %%%% b
x0=zeros(k,1);
x1=x0;
kmax=8000;
[iter1,MSE1,t1,SUDI1,E1,x1]=Algor2(R,f,b,x0,x1,kmax,k,p);
 E1,SUDI1(end), t1,
[iter2,MSE2,t2,E2,x2]=SCDKL(R,f,b,x0,x1,kmax,k,p);
 E2,t2,
[iter3,MSE3,t3,SUDI3,E3,x3]=SPK(R,f,b,x0,x1,kmax,k,p);
 E3,SUDI3(end),t3,
[iter4,MSE4,t4,SUDI4,E4,x4]=DTCRT(R,f,b,x1,kmax,k,p);
 E4,SUDI4(end),t4,
N=k;
figure(1)
scrsz = get(0,'ScreenSize');
set(1,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(5,1,1)
plot(f,'LineWidth',1.5)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Original signal (N= %g, number of nonzeros = %g)',N,p))
axis(v)

subplot(5,1,2)
plot(x1,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('Alg.2 (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter1,t1,MSE1))
axis(v)

subplot(5,1,3)
plot(x2,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('SCDKL (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter2,t2,MSE2))
axis(v)

subplot(5,1,4)
plot(x3,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('SPK (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter3,t3,MSE3))
axis(v)

subplot(5,1,5)
plot(x4,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('DTCRT (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter4,t4,MSE4))
axis(v)


