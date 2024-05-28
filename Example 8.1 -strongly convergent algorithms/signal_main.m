
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
lamda0=1;
[iter1,MSE1,t1,SUDI1,E1,x1]=Algor4(R,f,b,x0,x1,kmax,k,p);
E1, t1,
[iter2,MSE2,t2,SUDI2,E2,x2,F2]=suantai(R,f,b,x0,x1,kmax,k,p);
E2,t2,
[iter3,MSE3,t3,SUDI3,E3,x3,F3]=gibali(R,f,b,x0,x1,kmax,k,p);
E3,t3,
[iter4,MSE4,t4,SUDI4,E4,x4]=Lop(R,f,b,x1,kmax,k,p);
E4,t4,
[iter5,MSE5,t5,SUDI5,E5,x5]=SYGAAlg1(R,f,b,x0,x1,kmax,k,p,lamda0);
E5, t5,
[iter6,MSE6,t6,SUDI6,E6,x6]=DASAlg31(R,f,b,x0,x1,kmax,k,p);
E6, t6,
 N=k;
figure(1)
scrsz = get(0,'ScreenSize');
set(1,'Position',[10 scrsz(4)*0.1 0.9*scrsz(3)/2 3*scrsz(4)/4])
subplot(7,1,1)
plot(f,'LineWidth',1.5)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
set(gca,'FontName','Times')
set(gca,'FontSize',14)
title(sprintf('Original signal (N= %g, number of nonzeros = %g)',N,p))
axis(v)

subplot(7,1,2)
plot(x1,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('Alg.4(Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter1,t1,MSE1))
axis(v)

subplot(7,1,3)
plot(x2,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('SKCAlg.3.1 (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter2,t2,MSE2))
axis(v)


subplot(7,1,4)
plot(x3,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('GMVAlg.3.1 (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter3,t3,MSE3))
axis(v)

subplot(7,1,5)
plot(x4,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 N+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('LopezAlg.5.1 (Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter4,t4,MSE4))
axis(v)

subplot(7,1,6)
plot(x5,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 k+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('SYGAAlg.1(Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter5,t5,MSE5))
axis(v)

subplot(7,1,7)
plot(x6,'LineWidth',1.5)
set(gca,'FontName','Times')
set(gca,'FontSize',14)
top = max(f(:));
bottom = min(f(:));
v = [0 k+1 bottom-0.1*(top-bottom)  top+0.1*((top-bottom))];
title(sprintf('DASAlg.3.1(Iter=%5.3g, CPU=%5.3g, MSE=%5.3g)',iter6,t6,MSE6))
axis(v)






% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% µü´úÍ¼
% % linewidth = 2;
% % axesFontSize = 10;
% % labelFontSize = 10;
% % legendFontSize = 10;
% % 
% % resolution = 300; % output resolution
% % output_size = 300 *[10, 8]; % output size
% % figure(102), clf;
% % 
% % hold on,
% % p1f = semilogy(F1, 'b--', 'LineWidth',linewidth);
% % p2f = semilogy(F2, 'g-', 'LineWidth',linewidth);
% % p3f = semilogy(F3, 'r-', 'LineWidth',linewidth);
% % grid on;
% % ax.GridLineStyle = '--';
% % title({'\textbf{n=8000}'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% % ylb = ylabel({'$E_{n}$'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% % set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% % xlb = xlabel({'\vspace{-1.0mm}';'\textbf{Number of Iterations($n$)}'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% % set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);
% % lg = legend([p1f,p2f,p3f], ...
% %    'Algorithm 1','Algorithm 3.1[17]','Algorithm 2[27]');
% % set(lg,'FontSize', legendFontSize);
% % set(lg, 'Interpreter', 'latex');
% 
% % linewidth = 2;
% % axesFontSize = 10;
% % labelFontSize = 10;
% % legendFontSize = 10;
% % 
% % resolution = 300; % output resolution
% % output_size = 300 *[10, 8]; % output size
% % figure(2), clf;
% % 
% % hold on,
% % p1f = semilogy(SUDI1, 'b--', 'LineWidth',linewidth);
% % p2f = semilogy(SUDI2, 'g-', 'LineWidth',linewidth);
% % p3f = semilogy(SUDI3, 'r-', 'LineWidth',linewidth);
% % grid on;
% % ax.GridLineStyle = '--';
% % title({'\textbf{n=8000}'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% % ylb = ylabel({'$E_{n}$'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% % set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% % xlb = xlabel({'\vspace{-1.0mm}';'\textbf{Number of Iterations($n$)}'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% % set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);
% % lg = legend([p1f,p2f,p3f], ...
% %    'Algorithm 1','Algorithm 3.1[17]','Algorithm 2[27]');
% % set(lg,'FontSize', legendFontSize);
% % set(lg, 'Interpreter', 'latex');



