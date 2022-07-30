clear
clc
close all
%% 
BETA=0.2;
ALPHA=1;

INFINITY_ETA=10;
NODES=20;
D=INFINITY_ETA/(NODES-1);

F=0:0:0; G=F; H=G;
ETA=F;

F0=0; G0=0; GN=1;
SHOOT=0.001;
ADD=0.1*SHOOT;

F(1)=F0; G(1)=G0;
H(1)=SHOOT;

MAXIT=10000;
ITER=1;
ACCERR=0.0001;
ERROR=100;

while ((ERROR>ACCERR) && (ITER<MAXIT))
for I=1:NODES-1
    ETA(I+1)=I*D;
    F(I+1)=F(I)+D*G(I);
    G(I+1)=G(I)+D*H(I);
    H(I+1)=H(I)+D*(-ALPHA*F(I)*H(I) - BETA*(1-G(I)^2));
end

ERROR=abs(G(NODES)-GN);
H(1)=H(1)+ADD;

ITER=ITER+1;

end

for I=1:NODES
    OUT_M(:,1)=ETA;
    OUT_M(:,2)=F;
    OUT_M(:,3)=G;
    OUT_M(:,4)=H;
end
%% 
figure
plot (ETA,F,'linewidth',1.8)
hold on
plot (ETA,G,'linewidth',1.8)
plot (ETA,H,'linewidth',1.8)

xlabel ('\eta','FontSize',14,'color','b')
ylabel ('Functions','FontSize',14,'color','b')
title (['\beta = ',num2str(BETA)],'FontSize',12,'color','R')

grid on
axis square

legend('F(\eta)',"F^'(\eta)","F^'^'(\eta)",'FontSize',10,'color','G')

%% 
V=1.48E-5;
C=5;

LX=1;
LY=1;
DX=LX/(NODES-1);
DY=LY/(NODES-1);

X=0:0:0; Y=X; KSI=X;
S=zeros(NODES,NODES);
U_VEL=S; V_VEL=S;

for I=1:NODES
    X(I)=(I-1)*DX;
    KSI(I)=sqrt((V*(2-BETA))/C)*X(I)^((1-BETA)/(2-BETA));
    for J=1:NODES
        Y(J)=(J-1)*DY;
        
        INPUT=(Y(J)/(sqrt((2-BETA)*V/C)))*X(I)^(-(1-BETA)/(2-BETA));
        S(I,J)=sqrt(C*(2-BETA)*V) * X(I)^(1/(2-BETA)) * F_ETA(INPUT);
    end
end
S(1,:)=S(2,:);

for I=2:NODES-1
    for J=2:NODES-1
        U_VEL(1,J)=(S(1,J+1)-S(1,J-1))/(2*DY);
        U_VEL(I,J)=(S(I,J+1)-S(I,J-1))/(2*DY);
        
        U_VEL(I,1)=(S(I,2)-S(I,1))/(DY);
        U_VEL(I,NODES)=(S(I,NODES)-S(I,NODES-1))/(DY);
        
        V_VEL(1,J)=-(S(2,J)-S(1,J))/(DX);
        V_VEL(NODES,J)=-(S(NODES,J)-S(NODES-1,J))/(DX);
        
        V_VEL(I,J)=-(S(I+1,J)-S(I-1,J))/(2*DX);
        V_VEL(I,1)=-(S(I+1,1)-S(I-1,1))/(2*DX);
    end
end


figure
hold on
quiver(Y,X,U_VEL,V_VEL,'linewidth',1.,'color','r')
contour(X,Y,S,12,'linewidth',1.,'color','B','linestyle','-')
legend ('Veloity Field','StreamLines','FontSize',10,'color','Y')
xlabel ('X','FontSize',14,'color','b')
ylabel ('Y','FontSize',14,'color','b')
title (['Flow Over a Wedge','\beta = ',num2str(BETA)],'FontSize',12,'color','R')
axis square
grid on


figure
hold on
contour(X,Y,U_VEL,14,'linewidth',1.2,'color','R','linestyle','-')
contour(X,Y,V_VEL,14,'linewidth',1.2,'color','G','linestyle','-')
legend ('uVelocity IsoLines','vVelocity IsoLines','FontSize',10,'color','Y')
xlabel ('X','FontSize',14,'color','b')
ylabel ('Y','FontSize',14,'color','b')
title (['Flow Over a Wedge','\beta = ',num2str(BETA)],'FontSize',12,'color','R')
axis square
grid on

%% 
FETA=OUT_M(:,2);
function OUTPUT=F_ETA(INPUT)
       p1 =    0.006314;
       p2 =     -0.0871;
       p3 =      0.4416;
       p4 =     0.02019;
       p5 =   -0.006609;
       OUTPUT = p1*INPUT^4 + p2*INPUT^3 + p3*INPUT^2 + p4*INPUT + p5;
end

