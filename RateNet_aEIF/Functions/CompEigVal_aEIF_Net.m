function Lambs = CompEigVal_aEIF_Net(ModPar,ConPar,r0,rx,NullArray,flag)
% This function computes the eigenvalues of an aEIF network 
%
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% r0 = vector of steady-state activity
% rx = vector of external rates
% NullArray = 3D array containing nullclines
%             first index = consecutive values for activity
%             second index = populations
%             third index = index of population that is held constant
% flag = 0/1: compute eigenvalues (analytically/numerically),
%        2: compute slopes of nullclines ...
% Lambs = max. eigenvalues OR flag that denotes stability (-1: stable,
%         1: unstable)

if flag==0
    
    tm = ModPar(:,1)/1000.0;

    eps = ConPar.R;
    N = ConPar.NN;
    J = ConPar.J; 
    NTypes = length(N); % plus external excitatory pool
    df = zeros(NTypes-1,NTypes-1);

    for i=1:NTypes-1
        
        % influence of external pool
        mu = eps(i,NTypes)*N(NTypes)*J(i,NTypes)*rx(i)*tm(i);
        sig2 = J(i,NTypes)*mu;

        for j=1:NTypes-1
            mu = mu + tm(i)*eps(i,j)*N(j)*r0(j)*J(i,j);
            sig2 = sig2 + tm(i)*eps(i,j)*N(j)*r0(j)*J(i,j).^2;
        end 

        % compute Jacobian matrix
        CModPar = num2cell(ModPar(i,:));
        [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
        lb = min(EL,Vr) - 20.0;
        
        S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r0(i)/gL).*x + sf^2*exp((x-VT)./sf));
        
        for j=1:NTypes-1
            dm = tm(i)*eps(i,j)*N(j)*J(i,j);
            ds2 = dm*J(i,j);
            if i==j
                dw = b*(tw/1000)/gL;
                q = integral2(@(V,u) exp(S(V)-S(u)).*(0.5*ds2*(S(u)-S(V)) - (dm-dw).*(u-V)),lb,Vup,@(V) max(V,Vr),Vup);
            else
                q = integral2(@(V,u) exp(S(V)-S(u)).*(0.5*ds2*(S(u)-S(V)) - dm.*(u-V)),lb,Vup,@(V) max(V,Vr),Vup);
            end
            df(i,j) = (ds2/sig2) * r0(i) - 4*(tau/1000.0)*r0(i)^2/sig2^2 * q;
        end   
    end
    
    df = df - eye(NTypes-1);
    df = diag(1./tm)*df;
    e = eig(df);
    Lambs = max(real(e));
    
elseif flag==1
    
    tm = ModPar(:,1)/1000.0;

    eps = ConPar.R;
    N = ConPar.NN;
    J = ConPar.J; 
    NTypes = length(N); % plus external excitatory pool
    df = zeros(NTypes-1,NTypes-1);

     for i=1:NTypes-1
         
         for j=1:NTypes-1
             
            r = r0;
            r_test = linspace(r(j)-0.05,r(j)+0.05,5);
            F_test = zeros(1,length(r_test));
            
            for k = 1:length(r_test)
                
                r(j) = r_test(k);
         
                mu = eps(i,NTypes)*N(NTypes)*J(i,NTypes)*rx(i)*tm(i);
                sig2 = J(i,NTypes)*mu;

                for l=1:NTypes-1
                    mu = mu + tm(i)*eps(i,l)*N(l)*r(l)*J(i,l);
                    sig2 = sig2 + tm(i)*eps(i,l)*N(l)*r(l)*J(i,l).^2;
                end 
                
                % compute transfer function  
                CModPar = num2cell(ModPar(i,:));
                [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
                lb = min(EL,Vr) - 20.0;

                S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(i)/gL).*x + sf^2*exp((x-VT)./sf));
                q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
                F_test(k) = 1000.0*sig2/(2*tau*q);
               
            end
            
            p = polyfit(r_test,F_test,1);
            df(i,j) = p(1);
         end
     end
    
     df = df - eye(NTypes-1);
     df = diag(1./tm)*df;
     e = eig(df);
     Lambs = max(real(e));       
       
else
    
    x1 = NullArray(:,1,1);
    y1 = NullArray(:,2,1);
    x2 = NullArray(:,1,2);
    y2 = NullArray(:,2,2);
    [~,id1] = min((x1-r0(1)).^2);
    [~,id2] = min((x2-r0(1)).^2);
    
    m1 = polyfit(x1(id1-1:id1+1),y1(id1-1:id1+1),1);
    m2 = polyfit(x2(id2-1:id2+1),y2(id2-1:id2+1),1);
    
    if m1(1)<m2(1)
        Lambs = -1;
    else
        Lambs = 1;
    end
    
end