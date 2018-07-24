function F = Zustandsgleichungen(x2,u,Flag)
% Die Formulierung der Zustandsgleichungen und ihre Ableitungen

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

% Hilsfunktionen
phi1 = a1*u.^2+b1*u+c1;
phi3 = a3*u.^2+b3*u+c3;
phiw = 10^(-8);
fw = 1./(x2-1.001*(1-a)).^3 + 1./(x2+1.001*(1-a)).^3;
Dfw = -3./(x2-1.001*(1-a)).^4 - 3./(x2+1.001*(1-a)).^4;
tau = s+t*(a./(1.001-a-x2));
ksi = 3*pi*Re*a * tau;  

%Dx3 = 0.5*u.^2;
if Flag == 0 
    %State equations
    Dx1 = -x2.^2+1 + (1./ksi).*u;
    Dx2 = 1./ksi.*(phi1.*x2 + phi3.*(x2.^3) + phiw.*fw);
    F = [Dx1;Dx2];
elseif Flag == 1
    %State equations
    Dx2 = 1./ksi.*(phi1.*x2 + phi3.*(x2.^3) + phiw.*fw);
    F = Dx2;
elseif Flag == 2
    ksix = (1/(3*pi*Re*a))*(-a*t./(a*t-s*(a+x2-1.001)).^2);
    DDx1 = -2*x2 + ksix .*u;
     Dx2_2 = (phi1.*x2 + phi3.*(x2.^3) + phiw.*fw);
    DDx2 = 1./ksi .*(phi1 + 3 *phi3.*(x2.^2) + phiw.*Dfw)+ksix .* Dx2_2;
    F = [DDx1;DDx2];
elseif Flag == 3
    Dx1 = -x2.^2+1 + (1./ksi).*u;
    Dx2 = 1./ksi.*(phi1.*x2 + phi3.*(x2.^3) + phiw.*fw);
    DDu = 1./ksi .* ((2*a1*u + b1).*x2 + (2*a3*u+b3).*(x2.^3));
    %F = [DDu Dx1;1./ksi Dx2];
    F = [DDu;ksi];
end