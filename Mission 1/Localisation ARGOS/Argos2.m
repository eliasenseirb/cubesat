clear all
close all
clc 

%% Argos 2ème version : 

RE = 6378.137e3;                                  % Taille du demi grand axe en m
f = 1/298.257223563;                              % Aplatissement de l'ellipsoide
RP = RE*(1-f);  
Vs=7e3;                                           % Vitesse du satellite en m/s
c=physconst('LightSpeed');                        % Célérité de la lumière en m/s
hs=1500e3;                                         % Altitude du satellite en basse orbite
GE=RE+hs; 
GP=RP+hs;

ftk0 = 868e6;                   %Fréquence d'emission par la plateforme
phik0 = 44.8; %latitude initiale 
lambdak0 = -0.59; %longitude initiale
h=0; %Balise en mer par exemple
variation=5; %pourcentage maximum de variation entre chaque mesure




mk= 4;                                         % Nombre de mesures de fréquences (doit etre >=3 pour pouvoir avoir assez d'equations)
% on va prendre un echantillon de 4 mesures sur un passage 
xk=zeros(4,mk); %ensemble des mesures sur un passage.
xk(1,:)= [lambdak0 phik0 h ftk0];

% for i=2:4
%     var= (-variation + (2*variation) * rand(1))/100;
%     xk(i,:)=[lambdak0+var*lambdak0  phik0+var*phik0 h ftk0+var*ftk0];
% end
% 
% Jgk=Jacobien_H(lambdak0,phik0,h,ftk0); %Jacobienne pour la premiere mesure;
% 
% for i=2:4
%     Jgktmp=Jacobien_H(xk(i,1),xk(i,2),xk(i,3),xk(i,4));
%     Jgk=[Jgk; Jgktmp];
% end

% fr0= H(lambdak0,phik0,h,ftk0,-1);
% var= [ (-variation + (2*variation) * rand(1))/100  (-variation + (2*variation) * rand(1))/100  (-variation + (2*variation) * rand(1))/100]
% zk = [fr0 (fr0+var(1)*fr0) (fr0+var(2)*fr0) (fr0+var(3)*fr0)];     % Mesures de mk frequences reçues au k° passage satellite

sigma2k = 1;                                    % Variance du bruit
Rk = sigma2k*eye(mk);
