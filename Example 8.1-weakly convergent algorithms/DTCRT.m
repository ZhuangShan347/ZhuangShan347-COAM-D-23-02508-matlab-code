function [iter,MSE,t,SUDI,E,x1]=DTCRT(A,xs,b,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
%%%%==参数设置===
kappa=2e-4;
delta=0.499;
tau0=2;
%%%%开始迭代%%%%
while 1
    y=x1;
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
    r=tau0*norm(AA)-delta*norm(y-w);
    while r>0
        tau0=0.27*tau0;
        AA=A'*(A*y-A*w);
        r=tau0*norm(AA)-delta*norm(y-w);
    end
    d=y-w-tau0*AA;
    alp=(dot(y-w,d)+tau0*norm(A*w-b))/norm(d)^2;
    x1=y-kappa*alp*d;
    E=(norm(x1-xs))/(max(1,norm(x1)));
    MSE=(norm(x1-xs)^2)/k;
    if iter>kmax
        break;
    end
    iter=iter+1;
    SUDI(iter)=E;
end
t=toc;
end