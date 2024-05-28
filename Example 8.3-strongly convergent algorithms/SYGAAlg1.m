function [k,t,F,SUDI]=SYGAAlg1(A,x1,x0,lamda0)
clc
global epss
tic
k=1;
delta=0.3;
while 1
    x=x1-x0;
    if mod(k,2)==0
        y=x1;
    else
        y=x1+0.3*x;
    end
    x0=x1;
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
    grad_fn=A'*(y1-pq);%% 计算函数在 y_n处的梯度
    dd=lamda0* grad_fn;
    x=y-dd; %%% 接下来计算 C上的投影
    y2=y(2:10); %% n 在变化的时候，这里也在变化
    c=-y(1)+(y2'*y2);
    u=[-1; 2*y2];%%%grad_c
    u2=u'*u;
    c1=c-u'*dd; %% 出现了问题
    if c1>0
        x_bar=x-(c1*u)/u2;
    else
        x_bar=x;
    end
    y1=A*x_bar;   %%%AY
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
    grad_xx=A'*(y1-pq);%% 计算函数在 y_n处的梯度
    dist=lamda0*norm(grad_fn- grad_xx)-delta*norm(y- x_bar);
    while dist>0
        lamda0=0.5*lamda0;
        dd=lamda0* grad_fn;
        x=y-dd; %%% 接下来计算 C上的投影
        y2=y(2:10); %% n 在变化的时候，这里也在变化
        c=-y(1)+(y2'*y2);
        u=[-1; 2*y2];%%%grad_c
        u2=u'*u;
        c1=c-u'*dd; %% 出现了问题
        if c1>0
            x_bar=x-(c1*u)/u2;
        else
            x_bar=x;
        end
        y1=A*x_bar;   %%%AY
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
        grad_xx=A'*(y1-pq);%% 计算函数在 y_n处的梯度
        dist=lamda0*norm(grad_fn-grad_xx)-delta*norm(y- x_bar);
    end
    dd=lamda0*grad_xx;
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
    if  norm(x1-y)<epss
        break;
    end
    k=k+1;
    F(k)=norm(x1-x0);
    SUDI(k)=norm(x1-y);
end
t=toc;
end