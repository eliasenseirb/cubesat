function [cd_estime] = doppler_rate_esti(x,M,Np,chirp,T)
    %fonction permettant d'estimer le doppler rate
    
    temp=floor(length(x)/M); % Durée d'un chirp
    x=x(1:temp*M); % on redimensionne x pour le reshape
    
    Mat=reshape(x,[M,temp]); % on met en colonne les chirps
    for i=1:temp %calcul de R pour tous les chirps
        R(i) = 1/sqrt(M)*sum((Mat(:,i).*conj(chirp)').*exp(-1j*2*pi*(0:M-1)/M)');
        [maxi(i),ip(i)]=max(R(i));
    end
    ip=ip-1;
    for p=0:Np-2 %calcul de cdpl equation 2.51  
        for l=p+1:Np-1
            cd_pl(p+1,l) = 2/(T^2)*(ip(p+1+l)-ip(p+1))/l;  
        end
    end
    cd_estime = 0;
    
    for p=0:Np-2 %calcul equation 2.52
        somme_droite = 0;
        for l=p+1:Np-1
            somme_droite = somme_droite+cd_pl(p+1,l);
        end
        cd_estime=cd_estime+somme_droite/(Np-p+1);
    end
    cd_estime = cd_estime/(Np-1);
end

