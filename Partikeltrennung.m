function Partikeltrennung()
% Trennung zweier Partikel

%% Konstanten der Optimierung
% Anzahl der Zeitpunkte
N = 100;
% Anfangswerte der Steuerung
% Es wird über die Schrittweite und u gesteuert
u_0 = [zeros(1,N);ones(1,N)];
y0 = [0;0.2;0;0.2];

% Parameter für das Funktional
cu = 1/10;
cx = 20000;
cz = 1/100;

% Endwerte der Zustände
xT_1 = 0.3;
xT_2 = 0.6;
zT = 500;

% Beschränkungen an die Steuerung
lb = [-20*ones(1,N);zeros(1,N)];
ub = [20*ones(1,N);1000*ones(1,N)];
%% Zustandsgleichungen

a1 = -0.0052;
b1 = -0.0656;
c1 = 0.147;
a3 = 0.0071;
b3 = 0.0539;
c3 = -0.71;
% Reynoldzahl
Re = 10;
function erg = Zustandsgleichung(x,u,Flag,Radius)

if (Radius == 0.2)
    %Parameter des Modells
    a = 0.2;    
    s = 1.21;
    t = 0.126;
else 
    %Parameter des Modells
    a = 0.3;    
    s = 1.32;
    t = 0.192;
end
% Hilsfunktionen
phi1 =  a1*u.^2+b1*u+c1;
phi3 = a3*u.^2+b3*u+c3;
phiw = 10^(-8);
fw = 1./(x-1.001*(1-a)).^3 + 1./(x+1.001*(1-a)).^3;
Dfw = -3./(x-1.001*(1-a)).^4 - 3./(x+1.001*(1-a)).^4;
ksi = 3*pi*Re*a *(s+t*(a./(1.001-a-abs(x))));


if Flag == 'f'
    f1 = -x.^2 + 1 + (1./ksi).*u;
    f2 = 1./ksi.*(phi1.*x + phi3.*(x.^3) + phiw*fw);
    erg = [f1; f2];
elseif strcmp(Flag,'f1')
    f1 = -x.^2 + 1 + (1./ksi).*u;
    erg = f1;
elseif strcmp(Flag,'f2')
    f2 = 1./ksi.*(phi1.*x + phi3.*(x.^3) + phiw*fw);
    erg = f2;
elseif strcmp(Flag,'Df1')
    % Die Ableitung von 1/ksi nach x
    D1_ksi =(1/(3*pi*Re*a))*(-a*t./(a*t-s*(a+x-1.001)).^2);
    % Die Ableitung von f1 nach x
    Df1 = -2*x + D1_ksi.*u;
    erg = Df1;
elseif strcmp(Flag,'Df2')
    % Die Ableitung von 1/ksi nach x
    D1_ksi =(1/(3*pi*Re*a))*(-a*t./(a*t-s*(a+x-1.001)).^2);
    % Die Ableitung von f2 nach x
    f2 = 1./ksi.*(phi1.*x + phi3.*(x.^3) + phiw*fw);
    Df2 = 1./ksi .*(phi1 + 3*phi3.*(x.^2) + phiw*Dfw)+D1_ksi .* f2;
    erg = Df2;
end

end
%% Numerische Verfahren zur Lösung der DGL
% Expliziter Euler
function y = dgl(u,Flag)
    y = zeros(4,N+1);
    y(:,1) = y0;
    h = u(2,1);
        
    % Expliziter Euler
    if(Flag == 'e')
        for j=1:N
           y(1:2,j+1) = y(1:2,j) + h * Zustandsgleichung(y(2,j),u(1,j),'f',0.2);
           y(3:4,j+1) = y(3:4,j) + h * Zustandsgleichung(y(4,j),u(1,j),'f',0.3);
        end
    end
    if(Flag == 'i')
    % Impliziter Euler
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

% Das zu minimierende Zielfunktional
function J = Funktional(u)
   zx = dgl(u,'e');
   J = cx*(zx(4,end)-xT_2)^2+ cx*(zx(2,end)-xT_1)^2 +cu*u(2,1)*(u(1,:)*u(1,:)')+cz*(zx(3,end)- zT)^2+cz*(zx(1,end)- zT)^2;
end
 
options = optimoptions('fmincon','GradObj','off','Display','iter','Algorithm','sqp');
[u,fval,exitflag,output,lambda] = fmincon(@Funktional,u_0,[],[],[],[],lb,ub,@Zustandbeschr,options);

%% Zustandsbeschränkungen
% Beschränkung an die Zustände: x soll zwischen -0.8 und 0.8 liegen
function [c,ceq] = Zustandbeschr(u)
    ceq = [];
    y = dgl(u,'e');
    c = abs(y(2,:))-0.8;   
end

%% Plotten der optimalen Lösung
% Optimale Endzeit
T = sum(u(2,1))*N;
% Zustände mit optimaler Steuerung berechnen
y_opt = dgl(u,'e');

figure(2)
title('Expliziter Euler')

blub = u(2,1):u(2,1):T
figure(1)
plot(y_opt(1,:),y_opt(2,:),'g');
hold on
plot(y_opt(3,:),y_opt(4,:),'b');
xlabel('z')
ylabel('x')
legend('Partikel 1','Partikel 2')
figure(1)
plot(blub,u(1,:));
xlabel('t')
ylabel('u')


end
