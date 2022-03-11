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

% fr1 = ftk0 + 80;              %fréquence reçue au début d'un passage satellite
% fr2 = ftk0 + 80*15*60;        %fréquence reçue à la fin du même passage satellite

%phik0 = atan((GE^2*c*(fr1/ftk0 -1))/(sqrt(abs(Vs^2-c^2*(fr1/ftk0 -1)))));   % Expression de la longitude en fonction de fr1 en radian
phik0 = 44.8;

%lambdak0 = sin((GP^4*tan(phik0))/(GE^4)*sqrt(abs(Vs/c*(1-fr2/ftk0)-1)));   % Expression de la latitude en fonction de fr2 en radian
lambdak0 = -0.59; 
h=0; %Balise en mer par exemple

xk0=[lambdak0 phik0 h ftk0]; %valeur de départ 
variation=5;
    %Méthode Gauss-Newton (xk1_est = xk0_est + deltaxk0_est)

mk= 4;                                         % Nombre de mesures de fréquences (doit etre >=3 pour pouvoir avoir assez d'equations)
fr0= H(lambdak0,phik0,h,ftk0,-1);
var= [ (-variation + (2*variation) * rand(1))/100  (-variation + (2*variation) * rand(1))/100  (-variation + (2*variation) * rand(1))/100]
zk = [fr0 (fr0+var(1)*fr0) (fr0+var(2)*fr0) (fr0+var(3)*fr0)];     % Mesures de mk frequences reçues au k° passage satellite



sigma2k = 1;                                    % Variance du bruit
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

tabdxko=zeros(4,mk);
for i=1:mk
    J= Jacobien_H(lambdak0,phik0,h,ftk0)';

    if(zk(i)>=0)                                             % Effet doppler positif donc le satellite se rapproche de la balise
        gk0 = H(lambdak0,phik0,h,ftk0,1);
    else
        gk0 = H(lambdak0,phik0,h,ftk0,-1);
    end
    
    dxk0=inv(J'*inv(Rk)*J)*J'*inv(Rk)*(zk(i)-gk0);              % calcul de la petite variation pour raffiner l'estimation des coord
    tabdxko(:,i)=dxk0;
    xk1=xk0+dxk0;
    Xk_MAT(i+1,:)=xk1;
    lambdak0 = xk1(1);
    phik0 = xk1(2);
    ftk0 =xk1(4);
    xk0= [lambdak0 phik0 xk1(3) ftk0];
end

%Xk_MAT(:,1:2)=Xk_MAT(:,1:2)*180/pi;






