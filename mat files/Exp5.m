clc;
    zdata=[1 1 0 0 0.25
           2 1 2 0 0.1
           3 2 0 0 0.25
           2 1 3 0 0.1
           4 2 3 0 0.1];
ref=input('enter ref bus number');
cno=zdata(:,1);
fb=zdata(:,2);
tb=zdata(:,3);
r=zdata(:,4);
x=zdata(:,5);
n=max(max(fb),max(tb));
el=length(fb);
zbus=zeros(n,n);
for i=1:el
    p=fb(i);
    q=tb(i);
    z(i)=complex(r(i),x(i));
switch cno(i)
    case 1
        if(p~=ref)
            zbus(p,p)=z(i);
        else
            zbus(q,q)=z(i);
        end
disp('zbus after adding element');
disp(i); zbus
    case 2
        newb=max(p,q);
        oldb=min(p,q);
        for k=1:n
            zbus(k,newb)=zbus(k,oldb);
            zbus(newb,k)=zbus(oldb,k);
        end
        zbus(newb,newb)=zbus(newb,oldb)+z(i);
        disp('zbus after adding element');
disp(i); zbus
    case 3
        if (p==ref)
            oldb=q;
        end
        if (q==ref)
            oldb=p;
        end
        for k=1:n
            zb(k)=-zbus(oldb,k);
        end
        zbi=-zb(oldb)+z(i);
        for k=1:n
            for j=1:n
                zbus(k,j)=zbus(k,j)-zb(k)*zb(j)/zbi;
            end
        end
        disp('zbus after adding element'); disp(i);
        zbus
    case 4
        for k=1:n
            zb(k)=zbus(p,k)-zbus(q,k);
        end
        zbi=zb(p)-zb(q)+z(i);
        for k=1:n
            for j=1:n
                zbus(k,j)=zbus(k,j)-zb(k)*zb(j)/zbi;
            end
        end
        disp('zbus after adding element');
        disp(i);
        zbus
end
end
disp('zbus of the given network is:'); zbus