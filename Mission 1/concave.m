function [symb,Rw] = concave(z,Dp_estimes)
    a_tab = Dp_estimes -1;
    b_tab = Dp_estimes +1 ;
    
    p=10;
    N=2^p;
    
    for k=1:length(a_tab)
        points=linspace(a_tab(k),b_tab(k),N);
        for i=1:p
            mid1=a_tab(k)+(b_tab(k)-a_tab(k))/2;
            mid2=a_tab(k)+(b_tab(k)-a_tab(k))/2 + 1/N;
            R_w = [r(z(:,k),a_tab(k)),r(z(:,k),mid1),r(z(:,k),mid2),r(z(:,k),b_tab(k))];
            [~, indice] = max(R_w);
            if indice <=2
                b_tab(k)=mid1;
            else
                a_tab(k)=mid2;
            end
        end
        [Rw(k),sym(k)] = max([r(z(:,k),a_tab(k)),r(z(:,k),b_tab(k))]);
        if sym(k)==1
            symb(k)=a_tab(k);
        else
            symb(k)=b_tab(k);
        end
    end

end