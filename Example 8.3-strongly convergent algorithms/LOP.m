function [k,t,F,SUDI]=LOP(A,x1,x0)
clc
global epss
tic
k=1;
while 1
    y=x1;
    y1=A*y;   %%%AY
    yy=y1(2:20);
    Q=y1(1)+(yy'*yy)-1;%%%Q(AY)
    v=[1;2*yy]; %%% 计算出梯度
%     normv=norm(v);
    v2=v'*v;
    if Q>0
        pq=y1-Q*v/v2;
    else              %%%%%%%%% compute AY project Qn；
        pq=y1;
    end
    f=(y1-pq)'*(y1-pq);
    grad_fn=A'*(y1-pq);%% 计算函数在 y_n处的梯度
    ff=norm(grad_fn)^2;
    lamda=2*f/ff;
    dd=lamda* grad_fn;
    x=y-dd; %%% 接下来计算 C上的投影
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c-u'*dd; %% 出现了问题
    if c1>0
        x1=x-(c1*u)/u2;
    else
        x1=x;
    end
     alpha=1/k;
     x1=alpha*x0+(1-alpha)*x1;
     if  norm(x1-y)<epss
        break;
     end
      k=k+1; 
      F(k)=norm(x1-x0);
      SUDI(k)=norm(x1-y);
end 
t=toc;
end