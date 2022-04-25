clc
clear
close all
x_L = 0;
x_R = 1;
L = x_R - x_L;
N = 500; % number of grid in time
dx = L/(N - 1);
x = linspace(x_L, x_R, N);
lamda = 0.001; % lamda = dt/dx
dt = lamda * dx;
t = 0: dt: lamda;
M = numel(t);
u = zeros(numel(x), numel(t));
%% initial condition
for i = 1: N
    if x(i) < 1
        u(i, 1) = sin(4*pi.*x(i));
    end
end
%% U minus and U plus
u(1, :) = 0; % left boundary condition
u(N, :) = 0; % right boundary condition
uA = u; uLF = u; uLW = u; uMC = u; uUP  = u;
%% plot section
% Analytical
for j = 1: M
    for i = 1: N
        uA(i,j) = sin(4*pi.*(x(i)-0.5*t(j)));
        
    end
end
% Lax Friedrichs
for j = 1: M - 1
    for i = 2: N - 1
        uLF(i, j + 1) = uLF(i, j) - lamda * (f_LF(uLF(i + 1, j), uLF(i, j), lamda) - f_LF(uLF(i, j), uLF(i - 1, j), lamda));
    end
end
% Lax Wendroff
for j = 1: M - 1
    for i = 2: N - 1
        uLW(i, j + 1) = uLW(i, j) - lamda * (f_LW(uLW(i + 1, j), uLW(i, j), lamda) - f_LW(uLW(i, j), uLW(i - 1, j), lamda));
    end
end
% MacCormack
% step 1
for j = 1: M - 1
    for i = 2: N - 1
        uMC(i, j + 1) = uMC(i, j) - dt * ((1/(2*dx)*(sign(uMC(i,j))+1).*(f_MC1(uMC(i, j)) - f_MC1(uMC(i-1, j))))+...
            (1/(2*dx)*(1-sign(uMC(i,j))).*(f_MC1(uMC(i+1, j)) - f_MC1(uMC(i, j)))));
    end
end
% step 2
for j = 1: M - 1
    for i = 2: N - 1
        uMC(i, j + 1) = uMC(i, j) - lamda * (f_MC2(uMC(i + 1, j), uMC(i, j), lamda) - f_MC2(uMC(i, j), uMC(i - 1, j), lamda));
        
    end
end
% Upwind
for j = 1: M - 1
    for i = 2: N - 1
        uUP(i, j + 1) = uUP(i, j) - lamda * (f_UP(uUP(i, j)) - f_UP(uUP(i-1, j)));
        %uUP(i, j + 1) = uUP(i, j) - dt * ((1/(2*dx)*(sign(uUP(i,j))+1).*(f_UP(uUP(i, j)) - f_UP(uUP(i-1, j))))+...
            %(1/(2*dx)*(1-sign(uUP(i,j))).*(f_UP(uUP(i+1, j)) - f_UP(uUP(i, j)))));
    end
end
figure(1)
%plot(x, u(:, M), 'r--o', x, u(:, 1), 'b', 'LineWidth', 2, 'MarkerSize', 5)
plot(x, uA(:, M), 'r', x, uLF(:, M), 'b',x, uLW(:, M), 'g',x, uMC(:, M), 'm--'...
    ,x, uUP(:, M), 'y--', 'LineWidth', 2, 'MarkerSize', 5)
axis([0 1 -1.2 1.2])
grid
xlabel('x')
ylabel('u')
title('lamda = 0.001 , N = 500')
legend('A','LF','LW','MC','UP')

%% Fourier 
fA = fit(x',uA(:, M),'fourier8');
fLF = fit(x',uLF(:, M),'fourier8');
fLW = fit(x',uLW(:, M),'fourier8');
fMC = fit(x',uMC(:, M),'fourier8');
fUP = fit(x',uUP(:, M),'fourier8');

figure (2)
p1 = plot(fA,x,uA(:, M),'r');
hold on
p2 = plot(fLF,x,uLF(:, M));
p3 = plot(fLW,x,uLW(:, M));
p4 = plot(fMC,x,uMC(:, M));
p5 = plot(fUP,x,uUP(:, M));
set(p1,'color','r', 'LineWidth', 2, 'MarkerSize', 5);
set(p2,'color','b', 'LineWidth', 2, 'MarkerSize', 5);
set(p3,'color','g', 'LineWidth', 2, 'MarkerSize', 5);
set(p4,'color','m', 'LineWidth', 2, 'MarkerSize', 5);
set(p5,'color','y', 'LineWidth', 2, 'MarkerSize', 5);
grid
axis([0 1 -1.2 1.2])
xlabel('x')
ylabel('u')
title('Fourier b, lamda = 0.001 , N = 500')
%%  Phase difference 
iuA = find (abs(uA(:, M)-(max (uA(:, M)))) < 0.0001,2);
iuLF = find (abs(uLF(:, M)-(max (uLF(:, M)))) < 0.0001,1);
iuLW = find (abs(uLW(:, M)-(max (uLW(:, M)))) < 0.0001,1);
iuMC = find (abs(uMC(:, M)-(max (uMC(:, M)))) < 0.0001,1);
iuUP = find (abs(uUP(:, M)-(max (uUP(:, M)))) < 0.0001,1);

pdLF =  x(iuLF) - x(iuA);
pdLW =  x(iuLW) - x(iuA);
pdMC =  x(iuMC) - x(iuA);
pdUP =  x(iuUP) - x(iuA);
% Domain dissipation
ddLF = max (uA(:, M)) - max (uLF(:, M));
ddLW = max (uA(:, M)) - max (uLW(:, M));
ddMC = max (uA(:, M)) - max (uMC(:, M));
ddUP = max (uA(:, M)) - max (uUP(:, M));

q = {'Phase difference' 'Domain dissipation'};
data = [pdLF(2)     ddLF
pdLW(2)     ddLW
pdMC(2)    ddMC
pdUP(2)  ddUP];
z = figure(3);
zz = uitable(z,'data',data,'columnname',q);


%% Fourier function
%{
syms z  
sum = 0;  
y = sin(4*pi*z);   %function you want 
a0 = (1/pi)*int(y,z,-pi,pi); 
for n = 1:3 
        %finding the coefficients 
    an = (1/pi)*int(y*cos(n*z),z,-pi,pi); 
    bn = (1/pi)*int(y*sin(n*z),z,-pi,pi);    
    sum = sum+(an*cos(n*z)+bn*sin(n*z));  
end 
sum = (a0./2) + sum;
figure (2)
fplot(z,y,[-pi,pi]);
axis ([0 1 -1 1])

grid 
hold on  
fplot(z,sum,[-pi,pi]); 
%}