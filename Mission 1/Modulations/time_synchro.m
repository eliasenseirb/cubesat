function [n_estime] = time_synchro(K_estime,Np,Nsw,x,M,Fse,time_upsampled,B,Ts)
    a=K_estime-1/2-1;
    b=K_estime+1/2-1;
    I=b-a;
    max_errors = 1;
    N_it = log2((b-a)*M/max_errors);
    for i =1:N_it
        if H_w(Np,Nsw,a,x,M,Fse,time_upsampled,B,Ts) > H_w(Np,Nsw,b,x,M,Fse,time_upsampled,B,Ts)
            b=b-I/2;
            n_estime = a*M;
        else
            a=a+I/2;
            n_estime = b*M;
        end
        I=I/2;
    end
end

     
function [H] = H_w(Np,Nsw,w,x,M,chirp,time_upsampled,B,Ts) % Np : nombres de chirps up, Nsw : nombres de chirps de synchro, w : indice ou calculer la fonction
% x :signal
    nb_colonnes = Np+Nsw-1;
    if w<0
        %x=[exp(1j*2*pi.*time_upsampled.*fc(time_upsampled,0,B,Ts)),x];
        x=[zeros(1,M),x]; % on rajoute des 0 devant si le préambule se situe dans le premier paquet de M chirps
        w=1+w;
    end
    vect = zeros(1,nb_colonnes);
    for i =1:nb_colonnes
        %figure,plot(abs(fft(x((i+w-1)*M*Fse+1:(i+w)*M*Fse))))
        vect(i) = max(abs(fft(conj(chirp).*x((i+w-1)*M+1:(i+w)*M)))); % calcul des max des fft des paquets de M chirps
    end
    H = sum(vect); % on somme tous les maxs trouvés dans le préambule
end