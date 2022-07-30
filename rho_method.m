M_entry=1.2;
TtC=0.1;
MaxLoc=1;
AOA=15;
slip_tol=0.01;
chord=100;
Thickness=chord*TtC;
halfangle1=abs(atan(Thickness/MaxLoc/2));
halfangle2=abs(atan(Thickness/2/(chord-MaxLoc)));
P_entry=101325;
Asur1=sqrt(MaxLoc^2+(Thickness/2)^2);
Asur2=sqrt((chord-MaxLoc)^2+(Thickness/2)^2);
Asur3=Asur1;
Asur4=Asur2;
if halfangle1>AOA
    %Surface 1
    Flag1=3;
    Flag2=1;
    var1=M_entry;
    var2=halfangle1-AOA;
    [M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
    Psur1=P_entry*p2p1;
    Fy1=-Psur1*Asur1*cos(halfangle1);
    Fx1=Psur1*Asur1*sin(halfangle1);
    %Surface 2
    Flag3=2;
    var3=M2;
    var4=halfangle2+halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur2=Psur1/P3P4;
    Fy2=-Psur2*Asur2*cos(halfangle2);
    Fx2=-Psur2*Asur2*sin(halfangle2);
    Msur2=M4;
    %Surface 3
    Flag1=3;
    Flag2=1;
    var1=M_entry;
    var2=halfangle1+AOA;
    [M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
    Psur3=P_entry*p2p1;
    Fy3=Psur3*Asur3*cos(halfangle1);
    Fx3=Psur3*Asur3*sin(halfangle1);
    %surface 4
    Flag3=2;
    var3=M2;
    var4=halfangle2+halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur4=Psur3/P3P4;
    Fy4=Psur4*Asur4*cos(halfangle2);
    Fx4=-Psur4*Asur4*sin(halfangle2);
    Msur4=M4;
    
    
elseif halfangle1<AOA
    %Surface 1
    Flag3=2;
    var3=M_entry;
    var4=AOA-halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur1=P_entry/P3P4;
    Fy1=-Psur1*Asur1*cos(halfangle1);
    Fx1=Psur1*Asur1*sin(halfangle1);
    %Surface 2
    Flag3=2;
    var3=M4;
    var4=halfangle2+halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur2=Psur1/P3P4;
    Fy2=-Psur2*Asur2*cos(halfangle2);
    Fx2=-Psur2*Asur2*sin(halfangle2);
    Msur2=M4;
    %Surface 3
    Flag1=3;
    Flag2=1;
    var1=M_entry;
    var2=halfangle1+AOA;
    [M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
    Psur3=P_entry*p2p1;
    Fy3=Psur3*Asur3*cos(halfangle1);
    Fx3=Psur3*Asur3*sin(halfangle1);
    %surface 4
    Flag3=2;
    var3=M2;
    var4=halfangle2+halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur4=Psur3/P3P4;
    Fy4=Psur4*Asur4*cos(halfangle2);
    Fx4=-Psur4*Asur4*sin(halfangle2);
    Msur4=M4;
    
    
    
elseif halfangle1==AOA
    %Surface 1
    Psur1=P_entry;
    Fy1=-Psur1*Asur1;
    Fx1=0;
    %Surface 2
    Flag3=2;
    var3=M2;
    var4=halfangle2+halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur2=Psur1/P3P4;
    Fy2=-Psur2*Asur2*cos(halfangle2);
    Fx2=-Psur2*Asur2*sin(halfangle2)
    Msur2=M4;
    %Surface 3
    Flag1=3;
    Flag2=1;
    var1=M_entry;
    var2=halfangle1+AOA;
    [M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
    Psur3=P_entry*p2p1;
    Fy3=Psur3*Asur3*cos(halfangle1);
    Fx3=Psur3*Asur3*sin(halfangle1);
    %surface 4
    Flag3=2;
    var3=M2;
    var4=halfangle2+halfangle1;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    Psur4=Psur3/P3P4;
    Fy4=Psur4*Asur4*cos(halfangle2);
    Fx4=-Psur4*Asur4*sin(halfangle2);
    Msur4=M4;
    
    
end
Fx=Fx1+Fx2+Fx3+Fx4;
Fy=Fy1+Fy2+Fy3+Fy4;
L=Fy*cos(AOA)-Fx*sin(AOA);
D=Fx*cos(AOA)+Fy*sin(AOA);
cl=L/(0.7*P_entry*M_entry^2*chord);
cd=D/(0.7*P_entry*M_entry^2*chord);
slipangle=0;
step=pi/180;
Flag1=3;
Flag2=1;
var1=Msur2;
var2=AOA+halfangle2+slipangle;
[M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
P5=Psur2*p2p1;
Flag3=2;
var3=Msur4;
var4=AOA-halfangle2+slipangle;
[M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
P6=Psur4/P3P4;
while abs(P6-P5)>slip_tol
    if P6>P5
    slipangle=slipangle+step;
    Flag1=3;
    Flag2=1;
    var1=Msur2;
    var2=AOA+halfangle2+slipangle;
    [M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
    P5=Psur2*p2p1;
    
    Flag3=2;
    var3=Msur4;
    var4=AOA-halfangle2+slipangle;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    P6=Psur4/P3P4;
    if P5>P6
        step=step/10;
    end
    elseif P5>P6
    slipangle=slipangle-step;
    Flag1=3;
    Flag2=1;
    var1=Msur2;
    var2=AOA+halfangle2+slipangle;
    [M1,Beta,ro2ro1,p2p1,T2T1,M2,p02p01,Theta,Delta_S]=obliqueshock(Flag1,Flag2,var1,var2);
    P5=Psur2*p2p1;
    
    Flag3=2;
    var3=Msur4;
    var4=AOA-halfangle2+slipangle;
    [M3,M4,Mu3,Mu4,Theta2,T3T4,P3P4,ro3ro4]=fun3(var3,var4,Flag3);
    P6=Psur4/P3P4;
    if P6>P5
        step=step/10;
    end
            
    end
 
end
slipangle=slipangle;
wedge_xvector=[0 MaxLoc 100 MaxLoc 0];
wedge_yvector=[0 Thickness/2 0 -Thickness/2 0];
slip_xvector=[100 110 120 130];
slip_yvector=tan(slipangle).*slip_xvector-100*tan(slipangle);
chord_xvector=[-30 100];
chord_yvector=[0 0];
AOA_xvector=[-40 100];
AOA_yvector=tan(AOA).*AOA_xvector;

