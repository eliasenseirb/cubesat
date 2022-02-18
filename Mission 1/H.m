function [fr]=H(lambda,phi,h,ft,vs)
%fonction d'observation doppler de paramètres longitude lambda, altitude h,
%latitude phi et frequence de transmission du signal. 
Re = 6378.137e3 % Valeur du demi grand axe en km
f = 1/298.257223563 %  Valeur de l'applatissement
Rp = Re*(1-f); % Valeur du demi-petit axe 


vs=7e3; % Vitesse du satellite en m/s
GE=Re+h; 
GP=Rp+h;
theta = atan((GP/GE)*tan(phi)); % latitude paramétrique (ou réduite)

% Coordonnées cartésiennes de la balise à estimer.
x=GE*cos(theta)*cos(lambda); 
y=GE*cos(theta)*sin(lambda);
z=GP*sin(theta);

fr=ft*(1+1/c*(-vs*z/sqrt(x^2+y^2+z^2)));

end