function Optimierungfmincon()
% Es kann zwischen den Verfahren Expliziter und Impliziter Euler für
% die Lösung der Zustandsgleichung unterschieden werden
% Die Optimierung erfolgt mit fmincon SQP

% Wähle das Euler-Verfahren
% expliziter Euler: 'e'
% impliziter Euler: 'i'
C_Euler = 'e';

% Wähle ob mit oder ohne Zustandsbeschränkung
% Mit Zustandsbeschränkung: 
% - Penalty: 'P'
% - MATLAB: 'ZB'
% Ohne Zustandsbeschränkung: 'OB'
ZB = 'P';

% Wähle die Schranke für die Zustandsbeschränkung:
Grenze = 0.4;
%% Konstanten der Optimierung
% Anzahl der Zeitpunkte
N = 50;
% Anfangswerte der Steuerung
% Es wird über die Schrittweite und u gesteuert
u_0 = [zeros(1,N);ones(1,N)];

% Anfangswert des Zustandes
y0 = [0;0.3];

% Endwerte der Zustände
xT = 0.5;
zT = 200;

% Parameter für das Funktional
cu = 1/10;
cx = 20000;
cz = 1/100;

% Beschränkungen an die Steuerung
lb = [-20*ones(1,N);zeros(1,N)];
ub = [20*ones(1,N);1000*ones(1,N)];

%% Zustandsgleichungen
% Hier werden die Zustandsgleichungen formuliert und 
% Ableitungen der Zustandsgleichung berechnet.

%Parameter des Modells
a = 0.2;
a1 = -0.0052;
b1 = -0.0656;
c1 = 0.147;
a3 = 0.0071;
b3 = 0.0539;
c3 = -0.71;
s = 1.21;
t = 0.126;

% Reynoldzahl
Re = 10;

% Ruft die jeweilige Funktion oder Ableitung der Zustandsgleichung auf
function erg = Zustandsgleichung(x,u,Flag)

% Hilsfunktionen
phi1 = a1*u.^2+b1*u+c1;
phi3 = a3*u.^2+b3*u+c3;
phiw = 10^(-8);
fw = 1./(x-1.001*(1-a)).^3 + 1./(x+1.001*(1-a)).^3;
% Ableitung von fw
Dfw = -3./(x-1.001*(1-a)).^4 - 3./(x+1.001*(1-a)).^4;
ksi = 3*pi*Re*a *(s+t*(a./(1.001-a-abs(x))));

% Die rechte Seite
if Flag == 'f'
    f1 = -x.^2 + 1 + (1./ksi).*u;
    f2 = 1./ksi.*(phi1.*x + phi3.*(x.^3) + phiw*fw);
    erg = [f1; f2];
% Nur f^1 für den Zustand z    
elseif strcmp(Flag,'f1')
    f1 = -x.^2 + 1 + (1./ksi).*u;
    erg = f1;
% Nur f^2 für den Zustand x
elseif strcmp(Flag,'f2')
    f2 = 1./ksi.*(phi1.*x + phi3.*(x.^3) + phiw*fw);
    erg = f2;
% Ableitung von f^1
elseif strcmp(Flag,'Df1')
    % Die Ableitung von 1/ksi nach x
    D1_ksi =(1/(3*pi*Re*a))*(-a*t./(a*t-s*(a+x-1.001)).^2);
    % Die Ableitung von f1 nach x
    Df1 = -2*x + D1_ksi.*u;
    erg = Df1;
% Ableitung von f^2
elseif strcmp(Flag,'Df2')
    % Die Ableitung von 1/ksi nach x
    D1_ksi =(1/(3*pi*Re*a))*(-a*t./(a*t-s*(a+x-1.001)).^2);
    % Die Ableitung von f2 nach x
    f2_2 = (phi1.*x + phi3.*(x.^3) + phiw*fw);
    Df2 = 1./ksi .*(phi1 + 3*phi3.*(x.^2) + phiw*Dfw)+D1_ksi .* f2_2;
    erg = Df2;
end

end
%% Numerische Verfahren zur Lösung der DGL
% Expliziter Euler und Impliziter Euler

function y = dgl(u,Flag)
    y = zeros(2,N+1);
    y(:,1) = y0;
    h = u(2,1);
        
    % Expliziter Euler
    if(Flag == 'e')
        for j=1:N
           y(:,j+1) = y(:,j) + h * Zustandsgleichung(y(2,j),u(1,j),'f');
        end
    end
    % Impliziter Euler
    if(Flag == 'i')    
        for j=1:N
            x = y(2,j);
            max_iter = 100; tol = 1.0e-10;
            for k = 1 : max_iter       
                f = -x + y(2,j) + h * Zustandsgleichung(x,u(1,j),'f2');            
                Df = -1 + h*Zustandsgleichung(x,u(1,j),'Df2');
                if abs(f) < tol
                    %display('Nullstelle bestimmt');
                    break;
                elseif abs(Df) < tol
                    display('waagrechte Tangente'); break;
                end
                x = x - f/Df;
            end
            y(2,j+1) = x;
            y(1,j+1) = y(1,j) + h*Zustandsgleichung(x,u(1,j),'f1');
        end
    end  
end

%% Optimierung
function J = Funktional(u)
   % Zustände ermitteln
   zx = dgl(u,C_Euler);
   % Zielfunktional
   J = cx*(zx(2,end)-xT)^2 +cu*u(2,1)*(u(1,:)*u(1,:)')+cz*(zx(1,end)- zT)^2;
end

function J = Funktional_Penalty(u)
   % Zustände ermitteln
   zx = dgl(u,C_Euler);
   % Zielfunktional
   J = cx*(zx(2,end)-xT)^2 +cu*u(2,1)*(u(1,:)*u(1,:)')+cz*(zx(1,end)- zT)^2;
   % Falls wir das Penaltyverfahren verwenden
   J = 1000000*Penalty(u) + J;
end

options = optimoptions('fmincon','GradObj','off','Display','iter','Algorithm','sqp');
% Mit Zustandsbeschränkung
if(ZB == 'ZB')
    [u,fval,exitflag,output,lambda] = fmincon(@Funktional,u_0,[],[],[],[],lb,ub,@Zustandbeschr,options);
% Penalty-Verfahren
elseif(ZB == 'P')
        [u,fval,exitflag,output,lambda] = fmincon(@Funktional_Penalty,u_0,[],[],[],[],lb,ub,[],options);
% Ohne Zustandsbeschränkung
else
    [u,fval,exitflag,output,lambda] = fmincon(@Funktional,u_0,[],[],[],[],lb,ub,[],options);
    
end

    %% Zustandsbeschränkungen
% Beschränkung an die Zustände: x soll zwischen -0.4 und 0.4 liegen
function [c,ceq] = Zustandbeschr(u)
    ceq = [];
    y = dgl(u,C_Euler);
    c = abs(y(2,:))- Grenze;   
end
% Das Penalty-Verfahren
function g = Penalty(u)
    y = dgl(u,C_Euler);
    pos = abs(y(2,:));
    pos(pos<=Grenze) = 0;
    pos(pos>0) = pos(pos>0)-Grenze;
    g = sum(pos.^2);
end
%% Plotten der optimalen Lösung
% Optimale Endzeit
T = sum(u(2,1))*N;
% Zustände mit optimaler Steuerung berechnen
y_opt = dgl(u,C_Euler);

figure(2)

blub = u(2,1):u(2,1):T;
plot(y_opt(1,:),y_opt(2,:));
xlabel('z')
ylabel('x')
figure(1)
plot(blub,u(1,:));
xlabel('t')
ylabel('u')

%% Hamilton Funktion
% save('LagrangeOhneZB_EE')
% Adjungiertes System
p(:,N) =[2*cz*(y_opt(1,end)-zT);2*cx*(y_opt(2,end)-xT)];
for j = N-1:-1:1
    f_z = Zustandsgleichung(y_opt(2,j),u(1,j),'Df1');
    f_x = Zustandsgleichung(y_opt(2,j),u(1,j),'Df2');
    f_y = [0 f_z ; 0 f_x];
    p(:,j) = p(:,j+1) + h*(f_y'*p(:,j+1));
end
f = Zustandsgleichung(y_opt(2,1:N),u(1,:),'f');
H_p = cu*u(1,:)*u(1,:)' + p(1,:)*f(1,:)'+p(2,:)*f(2,:)'
end
