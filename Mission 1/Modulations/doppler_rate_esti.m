function [cd_estime] = doppler_rate_esti(x,M,Np,chirp,T)
    %fonction permettant d'estimer le doppler rate
    
    temp=floor(length(x)/M); % Durée d'un chirp
    x=x(1:temp*M); % on redimensionne x pour le reshape
    
    Mat=reshape(x,[M,temp]); % on met en colonne les chirps
    for i=1:temp
        R(i) = 1/sqrt(M)*sum(Mat(:,i).*chirp'.*exp(-1j*2*pi*(0:M-1)/M)');
        [maxi(i),ip(i)]=max(R(i));
    end
    
    for p=0:Np-2
        for l=p+1:Np-1
            cd_pl(p+1,l) = 2/(T^2)*(ip(p+1+l)-ip(p+1))/l;  
        end
    end
    cd_estime = 0;
    
    for p=0:Np-2
        somme_droite = 0;
        for l=p+1:Np-1
            somme_droite = somme_droite+cd_pl(p+1,l);
        end
        cd_estime=cd_estime+1/(Np-p)*somme_droite;
    end
    cd_estime = cd_estime/(Np-1);
end

