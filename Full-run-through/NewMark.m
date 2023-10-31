function u = NewMark(M, K, C, F, gamma, beta, Dt, N)
    %Newmark time integrator
    u=zeros(length(M),N);v=u;a=u;

    Mhat=M+C*gamma*Dt+K*Dt^2*beta;
    for n=1:N-1
        %predictor
        ap=a(:,n);
        vp=v(:,n)+Dt*((1-gamma)*a(:,n)+gamma*ap);
        up=u(:,n)+Dt*v(:,n)+Dt^2*((1/2-beta)*a(:,n)+beta*ap);
        %solution
        fhat=F(:,n+1)-(M*ap+C*vp+K*up);
        da=Mhat\fhat;
        %corrector
        a(:,n+1)=ap+da;
        v(:,n+1)=vp+Dt*gamma*da;
        u(:,n+1)=up+Dt^2*beta*da;
    end

end
