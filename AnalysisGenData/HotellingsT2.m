function [D,pv] = HotellingsT2(X,Y,alpha,flag)

% X - Matrix of measurments for first population: 
% each row represents an observation, each column represents a variable
% Y - Matrix of measurments for second population: 
% each row represents an observation, each column represents a variable
% alpha - significance level
% flag - 0 = two-sampled (independent), 1 = paired

    if flag==0
        Xm = mean(X);
        nX = size(X,1);
        Ym = mean(Y);
        nY = size(Y,1);
        p = size(X,2);

        SX = cov(X);
        SY = cov(Y);
        Sp = ((nX-1)*SX + (nY-1)*SY)/(nX+nY-2);

        T2 = (Xm-Ym)*inv(Sp)*(Xm-Ym)'/(1/nX + 1/nY);
        F = (nX+nY-p-1)*T2/(p*(nX+nY-2));
        %Ftable = finv(1-alpha,p,nX+nY-p-1);
        pv = 1-fcdf(F,p,nX+nY-p-1);

    elseif flag==1

        Z = X-Y;
        n = size(Z,1);
        p = size(Z,2);

        Zm = mean(Z);
        SZ = cov(Z);
        
        T2 = n*Zm*inv(SZ)*Zm';
        F = (n-p)*T2/(p*(n-1));
        %Ftable = finv(1-alpha,p,n-p);
        pv = 1-fcdf(F,p,n-p);
    end

    if pv < alpha
        D = 1;
    else
        D = 0;
    end

end
