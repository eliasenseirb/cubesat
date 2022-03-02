clc
clear;
close all;
 

%% Méthode des moindres carrés
    %Estimation des coordonnées initiales xk0=(lambdak0,phik0,ftk0)
ftk0 = 868e6;                   %Fréquence d'emission par la plateforme
phik0=44.8378;                  %Latitude de Bordeaux
lambdak0 = -0.594 ;             %Longitude de Bordeaux
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
zk = [5e3;6e3;7e3;8e3];     % Mesures de mk frequences à l'instant k

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









