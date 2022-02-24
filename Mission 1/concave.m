function [symb,Rw] = concave(z,Dp_estimes,M)
    a_tab = Dp_estimes -1;
    b_tab = Dp_estimes +1 ;
    
    p=12;
    N=2^p;
    
    for k=1:length(a_tab)
        for i=1:p
            mid1=a_tab(k)+(b_tab(k)-a_tab(k))/2;
            mid2=mid1 + 1/N;
            R_w = [r(z(:,k),a_tab(k),M),r(z(:,k),mid1,M),r(z(:,k),mid2,M),r(z(:,k),b_tab(k),M)];
            [~, indice] = max(R_w);
            if indice <=2
                b_tab(k)=mid1;
            else
                a_tab(k)=mid2;
            end
        end
        [Rw(k),sym(k)] = max([r(z(:,k),a_tab(k),M),r(z(:,k),b_tab(k),M)]);
        if sym(k)==1
            symb(k)=a_tab(k);
        else
            symb(k)=b_tab(k);
        end
    end

end