clc
clear;
close all;
 

%% Méthode des moindres carrés
    %Estimation des coordonnées initiales xk0=(lambdak0,phik0,ftk0)

RE = 6378.137e3;                                  % Taille du demi grand axe en m
f = 1/298.257223563;                              % Aplatissement de l'ellipsoide
RP = RE*(1-f);  
Vs=7e3;                                           % Vitesse du satellite en m/s
c=physconst('LightSpeed');                        % Célérité de la lumière en m/s
hs=1500e3;                                         % Altitude du satellite en basse orbite
GE=RE+hs; 
GP=RP+hs;

ftk0 = 868e6;                   %Fréquence d'emission par la plateforme

fr1 = ftk0 + 80;              %fréquence reçue au début d'un passage satellite
fr2 = ftk0 + 80*15*60;        %fréquence reçue à la fin du même passage satellite

phik0 = atan((GE^2*c*(fr1/ftk0 -1))/(sqrt(abs(Vs^2-c^2*(fr1/ftk0 -1)))));   % Expression de la longitude en fonction de fr1 en radian

disp(Vs^2-c^2*(fr1/ftk0 -1))
lambdak0 = sin((GP^4*tan(phik0))/(GE^4)*sqrt(abs(Vs/c*(1-fr2/ftk0)-1)));   % Expression de la latitude en fonction de fr2 en radian
disp(Vs/c*(1-fr2/ftk0)-1)

h=20;                           %Altitude moyen à Bordeaux
xk0=[lambdak0 phik0 h ftk0];

% deltafk0 = (rand()-0.5)*ftk0;         %Effet doppler initial (freq reçue - freq transmise)
% frk0 = deltafk0 + ftk0;
% 
% if deltafk0>=0
%     eloignement =1;
% else
%     eloignement=-1;
% end

    %Méthode Gauss-Newton (xk1_est = xk0_est + deltaxk0_est)
mk= 4;                      % Nombre de mesures de fréquences (doit etre >=3 pour pouvoir avoir assez d'equations)
zk = [ftk0 +80 ftk0 +72 ftk0+83 ftk0+81];     % Mesures de mk frequences reçues à l'instant k

sigma2k = 1;                % Variance du bruit
Rk = sigma2k*eye(mk);


% for j=1:mk
%     gk0(j,:)= H(lambdak0,phik0,h,H(lambdak0,phik0,h,ftk0,1)-zk(j),1);
% end

% J=zeros(mk,size(zk,1));
% for k=1:mk
%     J(k,:)=Jacobien_H(lambdak0,phik0,h,H(lambdak0,phik0,h,ftk0,1)-zk(k));
%     
% end

Xk_MAT = zeros(mk+1,4);
Xk_MAT(1,:) = xk0;

for i=1:mk
    J= Jacobien_H(lambdak0,phik0,h,ftk0)';
    gk0 = H(lambdak0,phik0,h,ftk0,1);
    
    dxk0=inv(J'*inv(Rk)*J)*J'*inv(Rk)*(zk(i)-gk0); % calcul de dxo, formule 1.11
    
    xk1=xk0+dxk0;
    Xk_MAT(i+1,:)=xk1;
    lambdak0 = xk1(1);
    phik0 = xk1(2);
    ftk0 =xk1(4);
end

Xk_MAT(:,1:2)=Xk_MAT(:,1:2)*180/pi;








