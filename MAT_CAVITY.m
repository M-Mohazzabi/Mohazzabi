clear
clc
close all
%% PARAMETERS 
H=1;

REYNOLDS=100;
U_0=1.0;
KIN_VISCOSITY=U_0*H/REYNOLDS;
PR=7;
K_F=KIN_VISCOSITY/PR;

N=40;
X=linspace(0,H,N);
Y=linspace(0,H,N);

DX=X(2)-X(1);
DY=Y(2)-Y(1);

B=DX/DY;

W=zeros(N,N);
S=W;
U=W;
V=W;
T=W;
T_STAR=T;
P=W;

RELAX=1.0;

U(:,1)=0;
U(:,N)=0;
U(1,:)=0;
U(N,:)=U_0;

V(:,1)=0;
V(:,N)=0;
V(1,:)=0;
V(N,:)=0;

T_COLD=10;
T_HOT=100;

T(:,1)=0;
T(:,N)=0;
T(1,:)=T_COLD;
T(N,:)=T_HOT;

%% SOLUTION 
TOL=1E-3;
ERR=TOL+1;
K=1;
KMAX=1000;
RESIDUAL=zeros(KMAX,7);

while ((ERR>TOL) && (K<KMAX))
    UOLD=U;
    VOLD=V;
    WOLD=W;
    SOLD=S;
    TOLD=T;
for I=2:N-1
    for J=2:N-1
        S(:,1)=S(:,2);
        S(:,N)=S(:,N-1);
        S(1,:)=S(2,:);
        S(N,:)=S(N-1,:) + U_0*DY;
        
        W(:,1)=-((S(:,3) - 2*S(:,2) + S(:,1))/(DY^2));
        W(:,N)=-((S(:,N) - 2*S(:,N-1) + S(:,N-2))/(DY^2));
        W(1,:)=-((S(3,:) - 2*S(2,:) + S(1,:))/(DX^2));
        W(N,:)=-((S(N,:) - 2*S(N-1,:) + S(N-2,:))/(DX^2));
        
        S(I,J) = (1-RELAX)*S(I,J) + RELAX*(W(I,J)*DX^2 + S(I+1,J) + ...
                           S(I-1,J) + S(I,J+1) + S(I,J-1))/(2*(1+B^2));
        U(I,J)=(S(I,J+1) - S(I,J-1))/(2*DY);
        V(I,J)=(S(I-1,J) - S(I+1,J))/(2*DX);
        
        LEFT=U(I,J)*((W(I+1,J) - W(I-1,J))/(2*DX)) + ...
             V(I,J)*((W(I,J+1) - W(I,J-1))/(2*DY));
        W(I,J)=(1-RELAX)*W(I,J) + RELAX*(((KIN_VISCOSITY/(DX^2))*...
               (W(I+1,J) + W(I-1,J) + (B^2)*(W(I,J+1) + W(I,J-1)))) - LEFT)/...
               ((2*KIN_VISCOSITY/(DX^2))*(1+B^2));
        
        T(:,1)=TOLD(:,2);
        T(:,N)=TOLD(:,N-1);
        T(I,J)=(1-RELAX)*T(I,J) + ...
                  RELAX*0.25*(T(I+1,J)+T(I-1,J)+T(I,J+1)+T(I,J-1));
    end
end
ERRU=(norm(U)-norm(UOLD))/norm(U);
ERRV=(norm(V)-norm(VOLD))/norm(V);
ERRW=(norm(W)-norm(WOLD))/norm(W);
ERRS=(norm(S)-norm(SOLD))/norm(S);
ERRT=(norm(T)-norm(TOLD))/norm(T);

ERR=max([ERRU,ERRV,ERRW,ERRS,ERRT]);

RESIDUAL(K,1:7)=[K,ERRU,ERRV,ERRW,ERRS,ERRT,ERR];
disp([K,ERR])
K=K+1;
end

RE_X=zeros(N,N);
NU_X=RE_X;
H_X=RE_X;
NUSELT_AVG=0:0:0;

for I=2:N-1
    for J=2:N-1
        P(I,1)=((W(I,2)-W(I,1))/(DY));
        P(1,J)=((W(2,J)-W(1,J))/(DX));
        P(I,J)=((W(I+1,J)-W(I-1,J))/(2*DX));
        P(N,J)=((W(N,J)-W(N-1,J))/(DX));
        P(I,1)=((W(I,N)-W(I,N-1))/(DY));
    end
end
for I=2:N-1
    for J=1:N
        H_X(1,J)=K_F*((T(2,J)-T(1,J))/(DY))/(T_HOT-T_COLD);
        H_X(I,J)=K_F*((T(I+1,J)-T(I-1,J))/(2*DY))/(T_HOT-T_COLD);
        H_X(N,J)=K_F*((T(N,J)-T(N-1,J))/(DY))/(T_HOT-T_COLD);
    end
end
for I=1:N
    INTEGRAL=0;
    for J=1:N
        T_STAR(I,J)=(T(I,J)-T_COLD)/(T_HOT-T_COLD);
        NU_X(I,J)=H_X(I,J)*H/K_F;
        RE_X(I,J)=U(I,J)*H/KIN_VISCOSITY;
        
        INTEGRAL=INTEGRAL + NU_X(I,J)*DX;
    end
    NUSELT_AVG(I)=INTEGRAL/H;
end




%% PLOT 
figure
semilogx (RESIDUAL(1:K-1,1),RESIDUAL(1:K-1,2),'linewidth',2)
hold on
semilogx (RESIDUAL(1:K-1,1),RESIDUAL(1:K-1,3),'linewidth',2)
semilogx (RESIDUAL(1:K-1,1),RESIDUAL(1:K-1,4),'linewidth',2)
semilogx (RESIDUAL(1:K-1,1),RESIDUAL(1:K-1,5),'linewidth',2)
semilogx (RESIDUAL(1:K-1,1),RESIDUAL(1:K-1,6),'linewidth',2)
semilogx (RESIDUAL(1:K-1,1),RESIDUAL(1:K-1,7),'linewidth',2)

set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('Iteration','FontSize',14,'color','B')
ylabel ('Residuals','FontSize',14,'color','B')

title (['Re = ',num2str(REYNOLDS)],'FontSize',16,'color','K')
grid on
legend ('Residual_{(U)}','Residual_{(V)}','Residual_{(\omega)}',...
        'Residual_{(\psi)}','Residual_{(T)}',...
        'Residual_{MAX}','FontSize',14,'color','Y')
%%
figure
contourf(X,Y,U,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['u Velocity',' (Re = ',num2str(REYNOLDS),')'],'FontSize',16,'color','K')
grid on
axis square
%%
figure
contourf(X,Y,V,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['v Velocity',' (Re = ',num2str(REYNOLDS),')'],'FontSize',16,'color','K')
grid on
axis square
%%
figure
contourf(X,Y,W,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['Vorticity (\omega)',' (Re = ',num2str(REYNOLDS),')'],...
        'FontSize',16,'color','K')
grid on
axis square
%%
figure
contourf(X,Y,S,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['Stream Function (\psi)',' (Re = ',num2str(REYNOLDS),')'],...
        'FontSize',16,'color','K')
grid on
axis square
%%
figure
contourf(X,Y,T_STAR,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['NON-Dimensional Temperature (T^*)',' (Re = ',num2str(REYNOLDS),')'],...
        'FontSize',16,'color','K')
grid on
axis square
%%
figure
contourf(X,Y,P,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['Pressure (P)',' (Re = ',num2str(REYNOLDS),')'],...
        'FontSize',16,'color','K')
grid on
axis square
%%
figure
contourf(X,Y,NU_X,100,'edgecolor','none')
colorbar
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Y [m]','FontSize',14,'color','B')

title (['Local Nuselt Number (Nu_x)',' (Re = ',num2str(REYNOLDS),')'],...
        'FontSize',16,'color','K')
grid on
axis square
%%
figure
plot(Y,NUSELT_AVG,'-S','linewidth',2)
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('Y [m]','FontSize',14,'color','B')
ylabel ('Nu_{avg}','FontSize',14,'color','B')

title (['Average Nuselt Number Along Center Line',...
        ' (Re = ',num2str(REYNOLDS),')'],'FontSize',16,'color','K')
grid on
axis square
%%
figure
hold on
plot(X,RE_X(1,:),'-S','linewidth',2)
plot(X,RE_X(N,:),'-S','linewidth',2)
set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('X [m]','FontSize',14,'color','B')
ylabel ('Re_{x}','FontSize',14,'color','B')

title (['Local Reynolds Number Over Upper & Lower Walls',...
        ' (Re = ',num2str(REYNOLDS),')'],'FontSize',16,'color','K')
grid on
axis square
legend ('Lower Wall','Upper Wall','color','C')
%%
figure
hold on
plot(X,W(1,:),'-S','linewidth',2)
plot(X,W(N,:),'-S','linewidth',2)
plot(Y,W(:,1),'-S','linewidth',2)
plot(Y,W(:,N),'-S','linewidth',2)

set (gca,'xcolor','M','ycolor','M')
set (gca,'FontName','Times')

xlabel ('[X=Y] [m]','FontSize',14,'color','B')
ylabel ('Vorticity','FontSize',14,'color','B')

title (['Vorticity Over All Walls',...
        ' (Re = ',num2str(REYNOLDS),')'],'FontSize',16,'color','K')
grid on
axis square
legend ('Lower Wall','Upper Wall','Left Wall','Right Wall','color','C')
