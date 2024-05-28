clc
clear
close all
global epss
epss=1e-5;
n=10;
m=20;
aa=rand(m,n);
A=aa;
x0=rand(n,1);
x1=rand(n,1);
lamda0=0.01;
[k1,t1,F1,SUDI1]=suantai(A,x1,x0);
[k2,t2,F2,SUDI2]=gibali(A,x1,x0);
[k3,t3,F3,SUDI3]=PC(A,x1,x0,lamda0);
[k4,t4,F4,SUDI4]=LOP(A,x1,x0);
lamda0=2;
[k5,t5,F5,SUDI5]=SYGAAlg1(A,x1,x0,lamda0);
[k6,t6,F6,SUDI6]=DASAlg31(A,x1,x0);
k1,t1,k2,t2,k3,t3,k4,t4,k5,t5,k6,t6,
SKC=SUDI1(end),
GMV=SUDI2(end),
Alg3=SUDI3(end),   %% Algorithm 3 in example 8.3
LOP=SUDI4(end),
SYGAAlg1=SUDI5(end),
DASAlg31=SUDI6(end),










% F1(end);F2(end);F3(end);
% F4(end);F5(end);F6(end);F7(end);
%%%%%%%%%%%%% CPU Í¼
% hh = parula;
% %% distance error ||x_{k}-x^\star||  
% linewidth = 2;
% labelFontSize = 10;
% legendFontSize = 10;
% figure(101), clf;
% 
% hold on,
% p1d = semilogy((t1/(k1-1))*(0:(k1-1)),F1, 'g-', 'LineWidth',linewidth);
% p2d = semilogy((t2/(k2-1))*(0:(k2-1)),F2, 'b--', 'LineWidth',linewidth);
% uistack(p1d, 'bottom');
% grid on;
% ax = gca;
% ax.GridLineStyle = '--';
% % title({'\textbf{Sahu et al.Algorithm 3.1}'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% ylb = ylabel({'$\|x_{n+1}-x_{n}\|$'}, 'FontSize', labelFontSize,...
%     'FontAngle', 'normal', 'Interpreter', 'latex');
% 
% set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% xlb = xlabel({'\vspace{-1.0mm}';'\textbf{CPU}'}, 'FontSize', labelFontSize,...
%     'FontAngle', 'normal', 'Interpreter', 'latex');
% set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);
% 
% lg = legend([p1d,p2d],...
%     'Algorithm 3.1[17]','Algorithm 1');
% set(lg,'FontSize', legendFontSize);
% set(lg, 'Interpreter', 'latex');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% µü´úÍ¼
% linewidth = 2;
% axesFontSize = 10;
% labelFontSize = 10;
% legendFontSize = 10;
% 
% resolution = 300; % output resolution
% output_size = 300 *[10, 8]; % output size
% figure(102), clf;
% 
% hold on,
% p1f = semilogy(F1, 'r-', 'LineWidth',linewidth);
% p2f = semilogy(F2, 'k-', 'LineWidth',linewidth);
% p3f = semilogy(F3, 'b--', 'LineWidth',linewidth);
% p4f = semilogy(F4, 'b-', 'LineWidth',linewidth);
% p5f = semilogy(F5, 'c-', 'LineWidth',linewidth);
% p6f = semilogy(F6, 'k--', 'LineWidth',linewidth);
% p7f = semilogy(F7, 'm--', 'LineWidth',linewidth);
% 
% grid on;
% ax.GridLineStyle = '--';
% % title({'\textbf{Sahu et al.Algorithm 3.1}'}, 'FontSize', labelFontSize,...
% %     'FontAngle', 'normal', 'Interpreter', 'latex');
% ylb = ylabel({'$\|x_{n+1}-x_{n}\|$'}, 'FontSize', labelFontSize,...
%     'FontAngle', 'normal', 'Interpreter', 'latex');
% set(ylb, 'Units', 'Normalized', 'Position', [-0.1, 0.5, 0]);
% xlb = xlabel({'\vspace{-1.0mm}';'\textbf{Number of Iterations($n$)}'}, 'FontSize', labelFontSize,...
%     'FontAngle', 'normal', 'Interpreter', 'latex');
% set(xlb, 'Units', 'Normalized', 'Position', [1/2, -0.075, 0]);
% lg = legend([p1f,p2f,p3f,p4f,p5f,p6f,p7f], ...
%   'GMVAlg.3.1','SKCAlg.3.1','Alg.1','mxjCQ','LOP','SYGAAlg1', 'DASAlg1');
% set(lg,'FontSize', legendFontSize);
% set(lg, 'Interpreter', 'latex');



