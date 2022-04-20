function [f_c] = fc(t,gammap,B,Ts)
% fonction calculant les fréquences des symboles
    f_c=zeros(size(t));
    for k=1:length(t)
        if t(k)>=-Ts/2 && t(k)<-Ts/2+gammap
            f_c(k)=(B/(2*Ts))*(t(k)-2*gammap)+B;
        elseif t(k)>=-Ts/2+gammap && t(k)<=Ts/2
            f_c(k)=(B/(2*Ts))*(t(k)-2*gammap);
        end
    end
end

