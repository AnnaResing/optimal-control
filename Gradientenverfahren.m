%% Das Gradientenverfahren wird mithilfe des adjungierten Zustandes und 
%% des reduzierten Gradienten implementiert
%% Zusätzlich wird die Funktion Zustandsgleichungen.m benötigt
function Gradientenverfahren()
clear all
% Anzahl der Zeitschritte
N = 100;
% Schrittweite
h = 1/N;

% Anfangswerte der Steuerung
u = [zeros(1,N);0:1:N-1];

% Konstanten des Zielfunktionals
cu = 1/10;
cx = 50000;
cz = 1/100;

% Anfangswerte und Endwerte der Zustände
xT = 0.1;
zT = 200;
x0 = 0.6;
z0 = 0;

% Dann ergeben sich die Zustände mit dem expliziten Euler Verfahren
function y = euler(u)
    y = zeros(2,N+1);
    y(:,1) = [z0 ; x0];
    % Berechne die Zustände anhand der Steuerung mit dem expliziten Euler
    for j=1:N
        y(:,j+1) = y(:,j) + h * u(2,end)*Zustandsgleichungen(y(2,j),u(1,j),0);
    end
end
%% Schrittweitensteuerung nach Armijo

function lambda = Schrittweite(u,lambda0,d)
    stop = 0;
    lambda = lambda0;
    zx = euler(u);
    zx_neu = euler(u+[lambda(1,1) * d(1,:);lambda(2,1)*d(2,:)]);
    F = cx*(zx(2,end)-xT)^2 +cu*h*(u(1,:)*u(1,:)')+cz*(zx(1,end)- zT)^2;
    F_neu = cx*(zx_neu(2,end)-xT)^2 +cu*h*(u(1,:)*u(1,:)')+cz*(zx_neu(1,end)- zT)^2;
    % Für das erste lambda
    while (F_neu > F + 0.1*lambda(1,1)*dot(RedGr(1,:),d(1,:)'))
        lambda(1,1) = lambda(1,1)/2;
        zx = euler(u);
        zx_neu = euler(u+[lambda(1,1) * d(1,:);lambda(2,1)*d(2,:)]);
        F = cx*(zx(2,end)-xT)^2 +cu*h*(u(1,:)*u(1,:)')+cz*(zx(1,end)- zT)^2;
        F_neu = cx*(zx_neu(2,end)-xT)^2 +cu*h*(u(1,:)*u(1,:)')+cz*(zx_neu(1,end)- zT)^2;
        
    end
    % 2. Lambda
    while (F_neu > F + 0.1*lambda(2,1)*dot(RedGr(2,:),d(2,:)'))
        lambda(2,1) = lambda(2,1)/2;
        zx = euler(u);
        zx_neu = euler(u+[lambda(1,1) * d(1,:);lambda(2,1)*d(2,:)]);
        F = cx*(zx(2,end)-xT)^2 +cu*h*(u(1,:)*u(1,:)')+cz*(zx(1,end)- zT)^2;
        F_neu = cx*(zx_neu(2,end)-xT)^2 +cu*h*(u(1,:)*u(1,:)')+cz*(zx_neu(1,end)- zT)^2;
        
    end
end
%% Optimierung
k = 0;
flag = 0;
% Beginn des Gradientenverfahrens
while (k< 50000 && flag ==0) 

y = euler(u);

% Adjungierte Gleichung lösen
p(:,N+1) =[2*cz*(y(1,end)-zT);2*cx*(y(2,end)-xT)];

for j = N:-1:1
    temp = Zustandsgleichungen(y(2,j),u(1,j),2);
    f_x = [0 temp(1); 0 temp(2)];
    p(:,j) = p(:,j+1) + h*(f_x'*p(:,j+1))*u(2,end);
end

% Reduzierten Gradient berechnen
f_u = Zustandsgleichungen(y(2,1:N),u(1,:),3);
fp = (1./f_u(2,1:N).*p(1,2:N+1) + f_u(1,1:N) .* p(2,2:N+1));
f_v = Zustandsgleichungen(y(2,1:N),u(1,:),0);
fvp = f_v(1,1:N).*p(1,2:N+1) + f_v(2,1:N) .* p(2,2:N+1);
fup_final = fp + 2*cu*u(1,:);
fvp_final = fvp + cu*u(1,:).^2;
RedGr = [fup_final;fvp_final];

% Abstiegsrichtung
d = -RedGr;

% Abbruchbedingung einbauen
norm2 = norm(RedGr(1,:))+norm(RedGr(2,:))
if norm2 < 0.5
    flag = 1;
end

% Aktualisiere die Steuerung
u_last = u;
% Schrittweitensteuerung
%lambda = Schrittweite(u,ones(2,1),d)
% Feste Schrittweite
lambda = [0.01;0.01];
% Neues u
u = u + [lambda(1,1) * d(1,:);lambda(2,1)*d(2,:)];
% Nächste Iteration
k = k+1
end
%% Ergebnisse plotten

y = euler(u_last);
Ende = h * u_last(2,end);
F = cx*(y(2,end)-xT)^2 +cu*h*u_last(2,end)*(u_last(1,:)*u_last(1,:)')+cz*(y(1,end)- zT)^2;
figure(2)
plot(y(1,:),y(2,:));
xlabel('z')
ylabel('x')
figure(3)
plot(1:h*u_last(2,end):u_last(2,end),u_last(1,:));
xlabel('t')
ylabel('u')
end

