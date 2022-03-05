%Exercise 1 Mohammad Javad Mohazzabi
clc
clear

%% properties

tic %Time

T_0 = 50;  % base temperature
T_inf = 20; % Air temperature
u_0 = T_0 - T_inf;
H = 20; % convective heat transfer coefficient of air (W/m^2.K)
K = 237; % pure Aluminium conductive heat transfer coefficient (W/m.K)
r = 0.05; % radius
P = 2*pi*r; % Perimeter of the fin
a = pi*r^2; % cross section
m = sqrt((H*P)/(K*a));

%% domain & mesh

N=100;
L=1;
h=L/(N-1);
x=linspace(0,1,N);
A=zeros(N,N);
b=zeros(N,1);
u=zeros(N,1);

%% Left boundary
i=1;
A(i,i)=1;
b(i,1)=u_0;

%% Right boundary

i=N;
b(N,1)=0;
A(i,i-1)=1;
A(i,i)=((-(m^2*h^2)/2))-1-(h*H/K);

%% Internal domain

for i=2:N-1
    A(i,i)=(-(m^2*h^2))-2;
    A(i,i+1)=1;
    A(i,i-1)=1;
    b(i,1)=0;
end

%% Calculations and Results

AS=sparse (A);
u=gmres(AS,b,[],1e-5,N);
T_Kelvin = u+T_inf+273.15;
format long

%% Plot

plot(x,T_Kelvin,'linewidth',2)
xlabel('x')
ylabel('T')
legend(strcat('h=', num2str(h)))
hold on
toc