num=[10]
den=[0.04 0.5 1 0]
sys=tf(num,den)
bode(sys)
margin(sys)
[Gm,Pm,Wcg,Wcp]=margin(sys)