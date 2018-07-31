function PMP()
% Überprüfen der Bedingungen des Pontrjaginschen Maximumprinzip
% Optimale Zustände und Steuerung laden
load('LagrangeMitZB.mat')
% Gradient von g
Dg = [2*cz*(ysol(1,end)-zT);2*cx*(ysol(2,end)-xT)];

% Die Lagrangemultiplikatoren von MATLAB
lambda = [lambda.eqnonlin(1:N+1) (lambda.eqnonlin(N+2:2*(N+1))) lambda.ineqnonlin]';

% Adjungiertes System
p(:,N+1) =[2*cz*(ysol(1,end)-zT);2*cx*(ysol(2,end)-xT)];
for j = N:-1:1
    f_x = Zustandsgleichungen(ysol(2,j),ysol(3,j),2);
    f_y = [0 f_x(1,1) ; 0 f_x(2,1)];
    p(:,j) = p(:,j+1) + ysol(4,1)*(f_y'*p(:,j+1));
end

f = Zustandsgleichungen(ysol(2,:),ysol(3,:),0);
w0 = cu*ysol(3,end)*ysol(3,end)'
w1= p(1,:)*f(1,:)'
w2 = p(2,end)*f(2,end)'
%H_p = cu*u*u' + p(1,:)*f(1,:)'+p(2,:)*f(2,:)'
u = ysol(3,:);
lambda = [lambda(:,2) lambda(:,2:end)];
% Hamiltonfunktion
H_p = cu*u*u' + p(1,:)*f(1,:)'+ p(2,:)*f(2,:)' 
% Transversalitätsbedingung
H_T = cu*u(end)*u(end)' + p(1,end)*f(1,end)'+p(2,end)*f(2,end)' 
g_T = (ysol(1,end)-zT)^2 + (ysol(2,end)-xT)^2
Transv = H_T + g_T

%% Plotten der Ergebnisse
figure(7)
t = 0:ysol(4,1):ysol(4,1)*N;
plot(t,lambda(1,:),'b',t,lambda(2,:),'r',t,lambda(3,:),'g')
legend('Multiplikator für z','Multiplikator für x','Multiplikator für h')
%title('Reale Lagrangemultiplikatoren')
ylabel('Multiplikator')
xlabel('t')
figure(3)
t = 0:ysol(4,1):ysol(4,1)*N;
plot(t,p(2,:),'b',t,lambda(2,:),'g')
%title('Reale Lagrangemultiplikatoren')
ylabel('Multiplikator')
text(470,0,'Explizit')
text(370,2500,'MATLAB')
figure(4)
bla = lambda(2,:)* (-1/ysol(4,1));
plot(t,bla,'b',t,p(2,:),'g')
%title('Veränderte Lagrangemultiplikatoren')
xlabel('t')
ylabel('Multiplikator')
text(80, -5000,'Schrittweite:')
text(170, -5000,int2str(ysol(4,1)))
% subplot(3,1,3)
% plot(1:N,bla,'b',1:N,p(2,:),'g')
figure(5)
%subplot(2,1,1)
plot(t,p(1,:),'g',t,lambda(1,:),'b');
%title('Reale Lagrangemultiplikatoren')
ylabel('Multiplikator')
legend('Explizit','MATLAB')
%subplot(2,1,2)
figure(6)
blub = lambda(1,:)*(-1/ysol(4,1));
plot(t,blub,'b',t,p(1,:),'g');
%title('Veränderter Lagrangemultiplikatoren')
ylabel('Multiplikator')
xlabel('t')

figure(1)
subplot(2,1,1)
plot(ysol(1,:),ysol(2,:));
hold on
plot(271.6966,0:0.02:0.8,'r--')
hold on
plot(371.5121,0:0.02:0.8,'r--')
text(275,0.3,'Grenzintervall')
title('Steuerung und Zustand des Partikels')
xlabel('Länge des Leiters')
ylabel('x(t)')
subplot(2,1,2)

plot(ysol(1,:),[ysol(3,2) ysol(3,2:end)]);
hold on
plot(271.6966,0:0.05:2,'r--')
hold on
plot(371.5121,0:0.05:2,'r--')
text(275,1,'Grenzintervall')
xlabel('Länge des Leiters')
ylabel('Steuerung')
end