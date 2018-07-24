function FDM()
% Alles wird diskretisiert und die Steuerung erfolgt über x,z,u und die
% Schrittweite
% Mit zusätzlichen Zustandconstraints
clear all
% Anzahl der Zeitpunkte
N = 200;
% Anfangswerte und Endwerte der Zustände
x0 = 0.6;
xT = 0.1;
zT = 200;

% Konstanten des Zielfunktionals
cx = 20000;
cz = 0.01;
cu = 1/2;

% Beschränkungen an die Steuerung
lb = [zeros(1,N+1);-1*ones(1,N+1);-20*ones(1,N+1);0.2*ones(1,N+1)];
ub = [1000*ones(1,N+1);1*ones(1,N+1);20*ones(1,N+1);100*ones(1,N+1)];

% Anfangswerte für die Steuerung
% 1. Zeile: z 2.Zeile: x 3. Zeile: u 4.Zeile: h
y_0 = [zeros(1,N+1);x0*ones(1,N+1);zeros(1,N+1);ones(1,N+1)];

%% Das nichtlineare Gleichungssystem als Constraint formulieren

% D- Matix für Finite Differenzen
D = spdiags(ones(N,1)*[-1 1],-1:0,N,N);

% Constraints der Zustandsgleichungen
function [c,ceq] = constr1(y)
    %c = [];
    c = y(2,:)-1;
    z = y(1,2:N+1);
    x = y(2,2:(N+1));
    control = y(3,2:N+1);
    Schrittweite = y(4,1);
    f = Zustandsgleichungen(x,control,0);
    fz = f(1,:)';
    fx = f(2,:)';
    fx(1) = fx(1) + 1/Schrittweite * x0;
    Fz = [y(1,1) ; 1/Schrittweite *D * z' - fz];
    Fx = [y(2,1)- x0 ; 1/Schrittweite * D * x' - fx];
    ceq = [Fz Fx ];
end
%% Optimierung
% Zu minimierendes Funktional
Objective =@(y) cu*y(4,1)*sum(y(3,:).^2)+cx*(y(2,end)- xT)^2 +cz*(y(1,end)-zT)^2;

options = optimoptions('fmincon','Display','iter','Algorithm','sqp','TolCon',1e-12,'MaxFunEvals',120000);
[ysol,fval,exitflag,output,lambda] = fmincon(Objective,y_0,[],[],[],[],lb,ub,@constr1,options);

%% Ergebnisse Plotten

% Plotten
h1 = ysol(4,1);
Th = h1:h1:(N+1)*h1;
save('LagrangeOhneZB')

figure(2)
plot(ysol(1,:),ysol(2,:));
xlabel('z')
ylabel('x')
figure(3)
plot(ysol(1,:),ysol(3,:));
xlabel('z')
ylabel('u')

end
