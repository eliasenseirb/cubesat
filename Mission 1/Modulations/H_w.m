function [H] = H_w(Np,Nsw,w,x,M)

    temp=floor(length(x)/M); % Durée d'un chirp
    x=x(1:temp*M); % on redimensionne x pour le reshape
    
    sig=reshape(x,[M,temp]); % on met en colonne les chirps
    Y=abs(fft(sig)); % calcul du module de la fft pour chaque colonnes 
    Y2=Y;
    Y=Y(:); % vectorisation du du module de la fft
    nb_colonnes = Np+Nsw-1;
    for i =1:nb_colonnes
        if i == 1
            vect(i)=max(Y(w*M:ceil(w)*M));
        elseif i==nb_colonnes
            vect(i)=max(Y(floor(w+Np+Nsw-1)*M:(w+Np+Nsw-1)*M));
        else
            vect(i)=max(Y(ceil(w+i-2)*M:ceil(w+i-1)*M));
        end
        
    end
    %H = sum(max(Y2(:,w:w+Np+Nsw-1)));
    H = sum(vect);
    %H = sum(max(Y(w*M:(w+Np+Nsw-1)*M))); % application de la formule H(w)
    test=1;
end