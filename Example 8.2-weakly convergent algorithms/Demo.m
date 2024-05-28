
addpath('./HNO','./Images','./solvers')
clear all; close all; clc;
%I = double(imread('cameraman.pgm'));
%I = double(imread('pepper.png'));%%ѡ��1
%I = double(imread('chart.tiff'));
%I = double(imread('house.png')); %%ѡ��2
%I = double(imread('satellite.png'));
I = double(imread('house.png')); %%ѡ��3
%I = double(imread('mxj.jpg'));
I = I/255;

[P,center] = psfGauss([9,9],4); %%9X9,��׼����4��P�൱��ģ������
B = imfilter(I,P,'symmetric'); %% ʵ�����ɵ�ģ��ͼ��
randn('seed',314); %% seed ��ʾ�������������һ���ģ� ���� 1x 314 ����
Bobs =B+(1e-3)*randn(size(B)); %%�������ӵ�ģ��ͼ����%%Ҳ���ǹ۲�ͼ��
figure(1) 
%subplot(1,1,1); imshow(I,[]); %ԭʼͼ��
%subplot(1,1,1); imshow(Bobs,[],'Border','tight');
% %%%%%%%===== Alg1====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.chi=0.333;              % ���Բ�����ѡ��;
para.tau0=0.03;              % ��ʼ����;
para.delta=0.49;             % �����Ĳ���;
para.kappa=1;                % ���x_n+1��ϵ��;
para.Eps=1e-10;              % ͶӰ����С����
para.scale = 1;
para.detail = 0;
out1 =Alg1(I,Bobs,P,center,para,stop);
fprintf('=====Alg1 =========\r\n');
out1.Time(end),
out1.SNR(end),

% %%%%%%%===== SCDKL====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.Eps=1e-10;              % ͶӰ����С����
para.scale = 1;
para.detail = 0;
out3=SCDKL(I,Bobs,P,center,para,stop);
fprintf('=====SCDKL =========\r\n');
out3.Time(end),
out3.SNR(end),

% %%%%%%%===== SPK====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.chi=0.09;               % ���Բ�����ѡ��;
para.tau0=2;                 % ��ʼ����;
para.delta=0.49;             % �����Ĳ���;
para.kappa=1;                % ���x_n+1��ϵ��;
para.Eps=1e-10;              % ͶӰ����С����
para.scale = 1;
para.detail = 0;
out2 =SPK(I,Bobs,P,center,para,stop);
fprintf('=====SPK =========\r\n');
out2.Time(end),
out2.SNR(end),

% % %%%%%%%===== DTCRT====%%%%%%%%%%%%
stop.MAX = 100;     stop.rule = 'EPS';     stop.eps = 1e-3;
para.chi=0.09;               % ���Բ�����ѡ��;
para.tau0=2;                 % ��ʼ����;
para.delta=0.49;             % �����Ĳ���;
para.kappa=1.9;              % ���x_n+1��ϵ��;
para.Eps=1e-10;              % ͶӰ����С����
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


