function [iter,MSE,t,SUDI,E,x1]=SPK(A,xs,b,x0,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
%%%%==��������===
theta=0.27;
kappa=2e-4;
delta=0.499;
tau0=2;
%%%%��ʼ����%%%%
while 1
    x=x1-x0;
    m=norm(x)^2;
    epsilon=1/iter^2;
    m=epsilon/m;
    if m==0
        chi=theta;
    else
        chi=min( theta,m);
    end
    y=x1+chi*x;
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
        w=z; %%C�ϵ�ͶӰ
    end
    AA=A'*(A*y-A*w);
    r=tau0*norm(AA)-delta*norm(y-w);
    while r>0
        tau0=0.3*tau0;
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