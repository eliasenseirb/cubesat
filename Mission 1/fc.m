function [f_c] = fc(t,gammap,B,Ts)
    if t>=-Ts/2 && t<-Ts/2+gammap
        f_c=(B/Ts)*(t-gammap)+B;
    elseif t>=-Ts/2+gammap && t<=Ts/2
        f_c=(B/Ts)*(t-gammap);
    else 
        f_c=0;
    end
end

