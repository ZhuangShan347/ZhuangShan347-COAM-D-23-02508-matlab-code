function [iter,MSE,t,SUDI,E,x1]=Lop(A,xs,b,x1,kmax,k,pp)
E=(norm(x1-xs))/(max(1,norm(x1)));
iter=1;
MSE=(norm(x1-xs))^2/k;
SUDI(iter)=E;
while 1
    y=x1;
    x0=x1;
    ab=A*y-b;
    f=(ab)'*(ab);
    grad_f=A'*ab; %% P_Q(y)=b,
    lamda=1.999*f/norm(grad_f)^2;
    dd=lamda*grad_f;
    x1=y-dd; %%接下来计算在C上的投影
    ksi=sign(y);
    ksi2=ksi'*ksi;
    c=norm(y,1)-pp-ksi'*dd;
    if c>0
        x1=x1-c*ksi/ksi2;
    end
    alpha=0.0000001/iter;
    x1=alpha*x0+(1-alpha)*x1;
    E=(norm(x1-xs))/(max(1,norm(x1)));
    MSE=(norm(x1-xs)^2)/k;
    if iter>kmax
        break;
    end
    iter=iter+1;
    SUDI(iter)=E;
end
t=toc;












