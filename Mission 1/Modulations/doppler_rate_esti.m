function [cd_estime] = doppler_rate_esti(x,M,Np,chirp,T)
    %fonction permettant d'estimer le doppler rate
    
    temp=floor(length(x)/(M)); % Durée d'un chirp
    x=x(1:temp*M); % on redimensionne x pour le reshape
    
    Mat=reshape(x,[M,temp]); % on met en colonne les chirps
%     for i=1:Np %calcul de R pour tous les chirps up du preambule
%         for j=0:M-1
%             R(j+1,i) = 1/sqrt(M)*sum((Mat(:,i).*chirp').*exp(-1j*2*pi*(0:M-1)*j/M).');
%         end
%         [~,ip1(i)]=max(abs(R(:,i)));
%     end

    for i=1:Np
        Mat_Detect(:,i)=Mat(:,i).*chirp'; % On multiplie chaque colonne par le chirp brut conjugué
    end
    [~,test_ip1]=max(abs(fft(Mat_Detect)));!
    [ip,~] = concave(Mat_Detect,test_ip1-1,M); % on utilise l'algo de concavité pour améliorer la valeur de positionnement des maxs
    %cd_estime = 2/T^2*(ip(end)-ip(1))/length(ip); % calcul de l'estimation du DR entre la première valeur de préambule et la dernière
    cd_estime=0;
    for p=1:Np-1 %calcul equation 2.52
        somme_droite = 0;
        for l=p+1:Np
               somme_droite = somme_droite+2/(T^2)*((ip(l)-ip(p)))/((l-p)*10);
        end
        cd_estime=cd_estime+somme_droite/(Np-p);
    end
    cd_estime = cd_estime/(Np-1);
end