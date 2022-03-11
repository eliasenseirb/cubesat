function [n_estime] = time_synchro(K_estime,Np,Nsw,x,M)
    a=K_estime-1/2;
    b=K_estime+1/2;
    I=b-a;
    max_errors = 1;
    N_it = log2((b-a)*M/max_errors);
    for i =1:N_it
        if H_w(Np,Nsw,a,x,M) > H_w(Np,Nsw,b,x,M)
            b=b-I/2;
            n_estime = a*M;
        else
            a=a+I/2;
            n_estime = b*M;
        end
        I=I/2;
    end
end

