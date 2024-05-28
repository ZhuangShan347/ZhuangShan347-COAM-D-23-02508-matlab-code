function [iter,MSE,t,SUDI,E,x1]=Algor4(A,xs,b,x0,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
%%%%==参数设置
tau0=0.1;
kappa=1.5;
delta=0.49;
%%%%开始迭代
while 1
    x=x1-x0;
    normx=norm(x);
    m1=(iter+1)^2*normx;
    if x==0
        theta=0.3;
    else
        theta=min(0.3,1/m1);
    end
    y=x1+theta*x;
    x0=x1;
    ab=A*y-b;
    grad_h=A'*ab; %% P_Q(y)=b,
    T=tau0*grad_h;
    z=y-T;
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp-(ksi'*T);
    if c>0
        w=z-c*ksi/ksi2;
    else
        w=z; %%C上的投影
    end
    AA=A'*(A*y-A*w);
    N=dot(AA,y-w);
    NN=delta*norm(y-w)^2/N;
    p=1/(1+iter)^1.1+1;
    if N>0
        tau=min(NN,p*tau0);
    else
        tau=p*tau0;
    end
    d=y-w-tau*AA;
    tt=A'*(A*w-b);
    alp=(dot(y-w,d)+0.5*tau*dot(tt,tt))/norm(d)^2;
    z=y-kappa*alp*d;
    f=0.001*x0;
    x1=0.00000001/(iter+1)*f+(1-0.00000001/(iter+1))*z;
    E=(norm(x1-xs))/(max(1,norm(x1)));
    MSE=(norm(x1-xs)^2)/k;
    if iter>kmax
        break;
    end
    tau0=tau;
    iter=iter+1;
    SUDI(iter)=E;
end
t=toc;
end