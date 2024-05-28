function out3 = SCDKL(I,Bobs,P,center,para,stop)
% I ........The original image;
% Bobs......The observed image which is blurred and noisy
% P ........PSF of the blurring operator
% center ...A vector of length 2 containing the center of the PSF
% W ........A function handle. For an image X, W(X)  is  an orthogonal
%                                 transform of the image X.
% WT .......A function handle. For an image X, WT(X) is the inverse
%       (with respect to the operator W) othogonal transform of the image X
% sigma ...Regularization parameter
% para......Parameters structure
% stop.MAX .......maximum number of iterations
% para.fig ....1 if the image is shown at each iteration, 0 otherwise (Default=1)
% pars.BC .....boundary conditions. 'reflexive'  (default)  or 'periodic
% OUTPUT

% X_out ..... Solution of the problem 
%               min{||A(X)-Bobs||^2+ sigma\|Wx\|_1
% fun_all:Array containing all function values obtained in the FISTA method

[m,n]=size(Bobs);        Pbig = padPSF(P,[m,n]); %直接调用P矩阵

flag = exist('para');
if(flag && isfield(para,'fig'))  % display the figure recovered by algorithm
    fig = para.fig;
else
    fig = 0;
end
if(flag && isfield(para,'method')) % choose the method 'LADM' or 'ISM'
    method = para.method;
else
    method = 'correction';
end
if(flag && isfield(para,'detail')) %display the details of iterations
    detail = para.detail;
else
    detail = 0;
end
if(flag && isfield(para,'scale'))  % Scaling parameter for objective fun.
    scale = para.scale;
else
    scale = 1;
end

Eps  = para.Eps;   clear opts.Eps;       % 投影的最小参数


%%====== Corresponding to Reflexive boundary condition =======
trans = @(X)dct2(X);   itrans = @(X)idct2(X); %2 维离散余弦变换    %2维逆变换
% computng the eigenvalues of the blurring matrix
e1 = zeros(m,n);       e1(1,1) = 1;
Sbig = dct2(dctshift(Pbig,center))./dct2(e1);%% 矩阵A的生成，线性算子A的生成 

%%============================================================
% computing the two dimensional transform of Bobs
Btrans = trans(Bobs); %%观察矩阵二维离散余弦变换

%The Lipschitz constant of the gradient of 0.5*||A(X)-Bobs||^2
%  L = 1.5*max(max(abs(Sbig).^2));

% initialization
x = Bobs;   In = norm(I,'fro');%%% F范数：F范数是把一个矩阵中每个元素的平方求和后开根号
x0=x;x1=x; 
if strcmp(stop.rule,'fixed')
    out3.SNR = zeros(1,stop.MAX);    out3.Time = zeros(1,stop.MAX);
end
time = cputime;
for iter = 1: stop.MAX  
    normx=norm(x1-x0,'fro');
    normx2=normx^2;
    m1=iter^2*max(normx2,normx);
    if x1-x0==0
        chi=0.3;
    else
        chi=min(0.3,1/m1);
    end
    v=x1+chi*(x1-x0);
    x0=x1;   
    D = Sbig.*trans(v)-Btrans;
    if  norm(D,'fro')<=Eps
        Pv=Sbig.*trans(v);
        
    else
        Pv=x+(Eps/ norm(D,'fro')).*D;
    end
    
    D=Sbig.*trans(v)-trans(Pv);
    Gv=scale*itrans( conj(Sbig).*D);
    Gv = real(Gv);     % Fv的函数值
    eta=max(1, norm(Gv,'fro'));
    tau0=0.5./eta;
    switch method
        case 'correction'
            
     x1=v-tau0*Gv; %%%计算出y 的值
            
    end
    if strcmp(stop.rule,'fixed')
        out3.SNR(iter) = 20*log10( In/norm(x1-I,'fro') );
        out3.SNRR(iter) = 20*log10( In/norm(x-I,'fro') );
        out3.Time(iter) = cputime - time;
        out3.image = x1;  out3.iter = iter;
    elseif strcmp(stop.rule,'EPS')
        out3.err = norm(x1-x,'fro')/norm(x,'fro');
        out3.SNR(iter) = 20*log10( In/norm(x1-I,'fro') );
        out3.SNRR(iter) = 20*log10( In/norm(x-I,'fro') );
        out3.Time(iter) = cputime - time;
        if out3.err <= stop.eps ||  iter >= stop.MAX
            out3.image = x1;   out3.iter = iter;
            return;
        end
        
    elseif strcmp(stop.rule,'SNR')
        out3.err = norm(x1-x,'fro')/norm(x,'fro');
        out3.SNR(iter) = 20*log10( In/norm(x1-I,'fro') );
        out3.SNRR(iter) = 20*log10( In/norm(x-I,'fro') );
        out3.Time(iter) = cputime - time;
        if out3.SNR(end) >= stop.eps ||  iter >= stop.MAX
            out3.image = x1;   out3.iter = iter;
            return;
        end
    else
        error('Please refer to ISM.m for the stopping rule!!')
    end
    if (detail)
        % Compute the l1 norm of the wavelet transform and the function value
        % and store it in the function values vector fun_all if exists.
        t = sum(sum(abs(W(x1))));
        out3.fun_val = norm(Sbig.*trans(x1)-Btrans,'fro')^2 + sigma*t;
        % printing the information of the current iteration
        fprintf('iter=%3d   fun_val=%4.5f   error=%2.3e   SNR=%2.2f\n',iter,out3.fun_val,out3.err,out3.SNR(iter));
        % Displaying the picture recovered by the algorithm;
        if (fig)
            figure(820);  imshow(x1,[]);
        end
    end
    
    % Replace the last iterates;
    
end