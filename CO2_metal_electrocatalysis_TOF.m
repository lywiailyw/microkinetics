x=[-1.5:0.01:0.5];y=[-1:0.01:1];[X,Y]=meshgrid(x,y);format long
k1=36000*exp(-19.88*(Y-0.1));k1r=36000*exp(19.88*(Y-0.1));
k2=36000*exp(-19.88*(0.05+X-Y));k2r=36000*exp(19.88*(0.05+X-Y));
k3=10^13*exp(-39.75*(-X-0.15));k3r=2.01*10^7*ones(size(X));
syms a b c;
for i=1:201;
    for j=1:201;
        m=k1(i,j);
        d=-k1r(i,j)-k2(i,j);
        f=k2r(i,j);
        g=k3r(i,j);
        h=k2(i,j);
        n=-k2r(i,j)-k3(i,j);
        f1=m*a+d*b+f*c;
        f2=g*a+h*b+n*c;
        f3=a+b+c-1;
        g1=subs(f1);
        g2=subs(f2);
        g3=subs(f3);
        z=solve(g1,g2,g3,a,b,c);
        M(i,j)=double(z.a);
        N(i,j)=double(z.b);
        L(i,j)=double(z.c);
    end
end
Z=log10(22.2*(k2.*N-k2r.*L));
contourf(X,Y,Z,'EdgeColor','none')