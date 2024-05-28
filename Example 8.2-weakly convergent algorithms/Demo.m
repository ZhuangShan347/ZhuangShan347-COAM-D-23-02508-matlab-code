
addpath('./HNO','./Images','./solvers')
clear all; close all; clc;
%I = double(imread('cameraman.pgm'));
%I = double(imread('pepper.png'));%%选定1
%I = double(imread('chart.tiff'));
%I = double(imread('house.png')); %%选定2
%I = double(imread('satellite.png'));
I = double(imread('house.png')); %%选定3
%I = double(imread('mxj.jpg'));
I = I/255;

[P,center] = psfGauss([9,9],4); %%9X9,标准差是4，P相当于模糊算子
B = imfilter(I,P,'symmetric'); %% 实际生成的模糊图像
randn('seed',314); %% seed 表示产生的随机数是一样的， 产生 1x 314 数列
Bobs =B+(1e-3)*randn(size(B)); %%把噪音加到模糊图像上%%也就是观察图像
figure(1) 
%subplot(1,1,1); imshow(I,[]); %原始图像
%subplot(1,1,1); imshow(Bobs,[],'Border','tight');
% %%%%%%%===== Alg1====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.chi=0.333;              % 惯性参数的选择;
para.tau0=0.03;              % 初始步长;
para.delta=0.49;             % 步长的参数;
para.kappa=1;                % 求解x_n+1的系数;
para.Eps=1e-10;              % 投影的最小参数
para.scale = 1;
para.detail = 0;
out1 =Alg1(I,Bobs,P,center,para,stop);
fprintf('=====Alg1 =========\r\n');
out1.Time(end),
out1.SNR(end),

% %%%%%%%===== SCDKL====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.Eps=1e-10;              % 投影的最小参数
para.scale = 1;
para.detail = 0;
out3=SCDKL(I,Bobs,P,center,para,stop);
fprintf('=====SCDKL =========\r\n');
out3.Time(end),
out3.SNR(end),

% %%%%%%%===== SPK====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.chi=0.09;               % 惯性参数的选择;
para.tau0=2;                 % 初始步长;
para.delta=0.49;             % 步长的参数;
para.kappa=1;                % 求解x_n+1的系数;
para.Eps=1e-10;              % 投影的最小参数
para.scale = 1;
para.detail = 0;
out2 =SPK(I,Bobs,P,center,para,stop);
fprintf('=====SPK =========\r\n');
out2.Time(end),
out2.SNR(end),

% % %%%%%%%===== DTCRT====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.chi=0.09;               % 惯性参数的选择;
para.tau0=2;                 % 初始步长;
para.delta=0.49;             % 步长的参数;
para.kappa=1.9;              % 求解x_n+1的系数;
para.Eps=1e-10;              % 投影的最小参数
para.scale = 1;
para.detail = 0;
out4 =DTCRT(I,Bobs,P,center,para,stop);
fprintf('=====DTCRT =========\r\n');
out4.Time(end),
out4.SNR(end),

figure(1) 
subplot(1,1,1); imshow(out1.image,[],'Border','tight');
figure(2) 
subplot(1,1,1); imshow(out3.image,[],'Border','tight');
figure(3) 
subplot(1,1,1); imshow(out2.image,[],'Border','tight');
figure(4) 
subplot(1,1,1); imshow(out4.image,[],'Border','tight');



% fprintf('=====LOPPRO =========\r\n');
% fprintf('==iter==Time====SNR======iter====Time====SNR=======\r\n');
% fprintf('%d & %2.2f  & %2.2f & %2.2f \r\n\n\n',out.iter,out.Time(end), out.SNRR,out.SNR(end));


