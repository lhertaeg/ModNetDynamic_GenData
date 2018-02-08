function Lambs = CompEigVal_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r0,rx,flag)
% This function computes the eigenvalues of an aEIF network (1 Exc, 1 IN)
% after learning (that is, several excitatory pools have been generated and 
% respond to specific stimuli) 
%
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% LrnPar = structure containing all parameters related to learning
% r0 = vector of steady-state activity
% rx = vector of external rates
% flag = 0/1: compute eigenvalues (analytically/numerically)
% Lambs = max. eigenvalues

tm = ModPar(:,1)/1000.0;

f = LrnPar.f;
p = LrnPar.p;
m = LrnPar.m; % more (g_+)
l = 1 - f*(m-1)/(1-f); % less (g_-)
SDJ = LrnPar.SDJ;  
L = (1 + SDJ^2);

eps = ConPar.R;
N = ConPar.NN;
J = ConPar.J; 
NTypes = length(N); % plus external excitatory pool
df = zeros(4,4);
KEE = tm(1)*eps(1,1)*N(1)*J(1,1);

if flag==0

    %%%%% excitatory pools after learning %%%%%
    for i=1:3 % exc pools (+, sel, 0)

        % influence of external pool
        mu = eps(1,NTypes)*N(NTypes)*J(1,NTypes)*rx(1)*tm(1);
        sig2 = J(1,NTypes)*mu*L;

        % influence of remaining populations (actually just IN)
        for j=2:NTypes-1
            mu = mu + tm(1)*eps(1,j)*N(j)*r0(2+j)*J(1,j);
            sig2 = sig2 + L*tm(1)*eps(1,j)*N(j)*r0(2+j)*J(1,j).^2;
        end

        % influence of exc pools
        if i==1 % +
            mu = mu + (((p-2)*l + m)*f*r0(1) + f*l*r0(2) + (1-p*f)*l*r0(3))*KEE;
            sig2 = sig2 + L*(((p-2)*l^2 + m^2)*f*r0(1) + f*l^2*r0(2) + (1-p*f)*l^2*r0(3))*KEE*J(1,1);
        elseif i==2 % sel
            mu = mu + ((p-1)*f*l*r0(1) + f*m*r0(2) + (1-p*f)*l*r0(3))*KEE;
            sig2 = sig2 + L*((p-1)*f*l^2*r0(1) + f*m^2*r0(2) + (1-p*f)*l^2*r0(3))*KEE*J(1,1);
        elseif i==3 % 0
            mu = mu + ((p-1)*f*r0(1) + f*r0(2) + (1-p*f)*r0(3))*KEE;
            sig2 = sig2 + L*((p-1)*f*r0(1) + f*r0(2) + (1-p*f)*r0(3))*KEE*J(1,1);
        end

        % compute Jacobian matrix
        CModPar = num2cell(ModPar(1,:));
        [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
        lb = min(EL,Vr) - 20.0;

        S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r0(i)/gL).*x + sf^2*exp((x-VT)./sf));         

        for j=1:NTypes+1 % 3x excitatory pools + remaining ones (IN)       
            if j<=3
                [dp,dr] = CompDeriv(i,j,p,l,m,f);
                dm = KEE*dp;
                ds2 = L*KEE*J(1,1)*dr;
            else
                dm = tm(1)*eps(1,2)*N(2)*J(1,2);
                ds2 = L*dm*J(1,2);
            end
            if i==j
                dw = b*(tw/1000)/gL;
            else
                dw = 0;
            end 
            q = integral2(@(V,u) exp(S(V)-S(u)).*(0.5*ds2*(S(u)-S(V)) - (dm-dw).*(u-V)),lb,Vup,@(V) max(V,Vr),Vup);      
            df(i,j) = (ds2/sig2)*r0(i) - 4*(tau/1000.0)*r0(i)^2/sig2^2 * q; 
        end   

    end

    %%%%% interneuron population after learning %%%%%

    % influence of external pool
    mu = eps(2,NTypes)*N(NTypes)*J(2,NTypes)*rx(2)*tm(2);
    sig2 = J(2,NTypes)*mu*L;

    % influence of excitatory pools
    rE = f*(r0(2) + (p-1)*r0(1)) + (1-p*f)*r0(3);
    mu = mu + tm(2)*eps(2,1)*N(1)*rE*J(2,1);
    sig2 = sig2 + L*tm(2)*eps(2,1)*N(1)*rE*J(2,1).^2;

    % influence of IN-population
    mu = mu + tm(2)*eps(2,2)*N(2)*r0(4)*J(2,2);
    sig2 = sig2 + L*tm(2)*eps(2,2)*N(2)*r0(4)*J(2,2).^2;

    % compute Jacobian matrix
    CModPar = num2cell(ModPar(2,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;
    
    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r0(4)/gL).*x + sf^2*exp((x-VT)./sf));    
     
    for j=1:NTypes+1 % 3x excitatory pools + remaining ones (IN)
        dw = 0;
        KIE = tm(2)*eps(2,1)*N(1)*J(2,1);
        if j==1
            dm = f*(p-1)*KIE;
            ds2 = L*dm*J(2,1);
        elseif j==2
            dm = f*KIE;
            ds2 = L*dm*J(2,1);
        elseif j==3
            dm = (1-p*f)*KIE;
            ds2 = L*dm*J(2,1);
        else
            dm = tm(2)*eps(2,2)*N(2)*J(2,2);
            ds2 = L*dm*J(2,2);
            dw = b*(tw/1000)/gL;
        end
        q = integral2(@(V,u) exp(S(V)-S(u)).*(0.5*ds2*(S(u)-S(V)) - (dm-dw).*(u-V)),lb,Vup,@(V) max(V,Vr),Vup);      
        df(4,j) = (ds2/sig2)*r0(4) - 4*(tau/1000.0)*r0(4)^2/sig2^2 * q; 
    end   
    
else
    
    %%%%% excitatory pools after learning %%%%%
    for i=1:3 % exc pools (+, sel, 0)
       
        for j=1:NTypes+1
            
            r = r0;
            r_test = linspace(r0(j)-0.05,r0(j)+0.05,5);
            F_test = zeros(1,length(r_test));
            
            for k = 1:length(r_test)
                
                r(j) = r_test(k);

                % influence of external pool
                mu = eps(1,NTypes)*N(NTypes)*J(1,NTypes)*rx(1)*tm(1);
                sig2 = J(1,NTypes)*mu*L;

                % influence of remaining populations (actually just IN)
                mu = mu + tm(1)*eps(1,2)*N(2)*r(4)*J(1,2);
                sig2 = sig2 + L*tm(1)*eps(1,2)*N(2)*r(4)*J(1,2).^2;

                % influence of exc pools
                if i==1 % +
                    mu = mu + (((p-2)*l + m)*f*r(1) + f*l*r(2) + (1-p*f)*l*r(3))*KEE;
                    sig2 = sig2 + L*(((p-2)*l^2 + m^2)*f*r(1) + f*l^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);
                elseif i==2 % sel
                    mu = mu + ((p-1)*f*l*r(1) + f*m*r(2) + (1-p*f)*l*r(3))*KEE;
                    sig2 = sig2 + L*((p-1)*f*l^2*r(1) + f*m^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);
                elseif i==3 % 0
                    mu = mu + ((p-1)*f*r(1) + f*r(2) + (1-p*f)*r(3))*KEE;
                    sig2 = sig2 + L*((p-1)*f*r(1) + f*r(2) + (1-p*f)*r(3))*KEE*J(1,1);
                end
                
                % compute transfer function  
                CModPar = num2cell(ModPar(1,:));
                [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
                lb = min(EL,Vr) - 20.0;            

                S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(i)/gL).*x + sf^2*exp((x-VT)./sf));
                q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
                F_test(k) = 1000.0*sig2/(2*tau*q);
            end
            
            coef = polyfit(r_test,F_test,1);
            df(i,j) = coef(1);
        end
    end
    
    
    %%%%% interneuron population after learning %%%%%

    for j=1:NTypes+1
            
            r = r0;
            r_test = linspace(r0(j)-0.05,r0(j)+0.05,5);
            F_test = zeros(1,length(r_test));
            
            for k = 1:length(r_test)
                
                r(j) = r_test(k);
            
                % influence of external pool
                mu = eps(2,NTypes)*N(NTypes)*J(2,NTypes)*rx(2)*tm(2);
                sig2 = J(2,NTypes)*mu*L;

                % influence of excitatory pools
                rE = f*(r(2) + (p-1)*r(1)) + (1-p*f)*r(3);
                mu = mu + tm(2)*eps(2,1)*N(1)*rE*J(2,1);
                sig2 = sig2 + L*tm(2)*eps(2,1)*N(1)*rE*J(2,1).^2;

                % influence of IN-population
                mu = mu + tm(2)*eps(2,2)*N(2)*r(4)*J(2,2);
                sig2 = sig2 + L*tm(2)*eps(2,2)*N(2)*r(4)*J(2,2).^2;
                
                % compute transfer function  
                CModPar = num2cell(ModPar(2,:));
                [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
                lb = min(EL,Vr) - 20.0;            

                S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(4)/gL).*x + sf^2*exp((x-VT)./sf));
                q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
                F_test(k) = 1000.0*sig2/(2*tau*q);
            end
            
            coef = polyfit(r_test,F_test,1);
            df(4,j) = coef(1);
    end    
end

% Computing eigenvalues
tm = [tm(1),tm(1),tm(1),tm(2)];
df = df - eye(NTypes+1);
df = diag(1./tm)*df;
e = eig(df);
Lambs = max(real(e));
    

function [dp,dr] = CompDeriv(i,j,p,l,m,f)
        if (i==1 && j==1)
            dp = ((p-2)*l + m)*f;
            dr = ((p-2)*l^2 + m^2)*f;
        elseif (i==1 && j==2)
            dp = f*l;
            dr = f*l^2;
        elseif (i==1 && j==3)
            dp = (1-p*f)*l;
            dr = (1-p*f)*l^2;
        elseif (i==2 && j==1)
            dp = (p-1)*f*l;
            dr = (p-1)*f*l^2;
        elseif (i==2 && j==2)
            dp = f*m;
            dr = f*m^2;
        elseif (i==2 && j==3)
            dp = (1-p*f)*l;
            dr = (1-p*f)*l^2;    
        elseif (i==3 && j==1)
            dp = (p-1)*f;
            dr = (p-1)*f;
        elseif (i==3 && j==2)
            dp = f;
            dr = f;
        elseif (i==3 && j==3)
            dp = (1-p*f); 
            dr = (1-p*f);
        end