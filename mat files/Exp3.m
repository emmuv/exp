%YBUS formation USING INSPECTION METHOD without tap changing transformers
ldata = [ 1 1 2 0.02+0.06i 0.03i
2 1 3 0.08+0.24i 0.025i
3 2 3 0.06+0.18i 0.02i
4 2 4 0.06+0.18i 0.02i
5 2 5 0.04+0.12i 0.015i
6 3 4 0.01+0.03i 0.01i
7 4 5 0.08+0.24i 0.025i];
lno = ldata (:,1)
fb = ldata (:,2)
tb = ldata (:,3)
z = ldata (:,4)
ych = ldata(:,5)
ln =length(lno);
n=max (max(fb),max(tb));
y=1./z;
ybus = zeros(n,n);
for i =1:ln
p=fb(i);
q=tb(i);
ybus(p,p)=ybus(p,p)+y(i)+ych(i);
ybus(q,q)=ybus(q,q)+y(i)+ych(i);
ybus(p,q)=-y(i);
ybus(q,p)=-y(i);
end
disp ('YBUS OF GIVE SYSTEM IS :');
disp (ybus);