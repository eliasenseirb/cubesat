function [phi,lambda]=init_localisation(vs,fr,ft,hs,R,etha,alpha,phis)
% Fonction calculant l'intersection cone-sphere en cartésien et en
% géographique
%Entrées : vs :vitesse du satellite
%          fr :fréquence reçue au niveau du satellite
%          ft :fréquence transmise par la balise
%          hs : altitude du satellite 
%          R: rayon de la sphère (terre dans notre cas)
%          etha : angle entre l'axe joignant le sommet du cône et le centre
%          de la sphère avec l'axe Z
%          alpha : azimuth (angle) du Nord vers l'axe du cone
%          phis : latitude au point sous-sommet du cone

%Sorties: phi : latitude de la balise
%         lambda : longitude de la balise

%thèse utilisée : https://journals.ametsoc.org/view/journals/apme/10/3/1520-0450_1971_010_0607_tioaca_2_0_co_2.xml

c= physconst('Lightspeed'); %célérité de la lumière
a = c/vs*(fr/ft-1);
theta = atan(sqrt(1/a^2-1));

if(etha ==0)    %cas où le cone est vertical
    Z2 = (R+H)*sin(theta)^2 + sqrt(R^2- hs^2*sin(theta)^2)*cos(theta);
else
    Z2 = 


