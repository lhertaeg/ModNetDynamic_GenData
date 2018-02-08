function E = CompEnergy_aEIF_Net_Learn(rates,ModPar,ConPar,LrnPar,FP,rx,flg)
% This function computes the enery landscape of an aEIF network after 
% learning (along the selctive pop-activity)
%
% rates = vector of rates at which energy landscape is computed
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% LrnPar = structure containing all parameters related to learning
% FP = matrix of steady-state activity (1st row: spon, 2nd row: memo)
% rx = vector of external rates
% flg = 0: 1D projection along FP-vector, 1: equivalent 1D model (along
%       nullcline of remaining populations)
% E = "energy" landscape computed at values given in rates

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
dr = zeros(1,4);
fr = zeros(1,length(rates));
E = zeros(1,length(rates));
KEE = tm(1)*eps(1,1)*N(1)*J(1,1);

if flg==0
    
    qr = (rates(end)-FP(1,2))/(FP(2,2)-FP(1,2));   
    r2 = FP(1,:) + qr*(FP(2,:)-FP(1,:));

    for i=1:length(rates)

        qr = (rates(i)-FP(1,2))/(FP(2,2)-FP(1,2));   
        r = FP(1,:) + qr*(FP(2,:)-FP(1,:));

        for id=1:4

            % influence of external pool
            if id<4
                mu = eps(1,NTypes)*N(NTypes)*J(1,NTypes)*rx(1)*tm(1);
                sig2 = J(1,NTypes)*mu*L;    
            else
                mu = eps(2,NTypes)*N(NTypes)*J(2,NTypes)*rx(2)*tm(2);
                sig2 = J(2,NTypes)*mu*L;
            end

            % influence of inh pool
            if id<4
                mu = mu + tm(1)*eps(1,2)*N(2)*r(4)*J(1,2);
                sig2 = sig2 + L*tm(1)*eps(1,2)*N(2)*r(4)*J(1,2).^2;
            else
                mu = mu + tm(2)*eps(2,2)*N(2)*r(4)*J(2,2);
                sig2 = sig2 + L*tm(2)*eps(2,2)*N(2)*r(4)*J(2,2).^2;
            end

            % influence of exc pools
            if id==1 % +
                mu = mu + (((p-2)*l + m)*f*r(1) + f*l*r(2) + (1-p*f)*l*r(3))*KEE;
                sig2 = sig2 + L*(((p-2)*l^2 + m^2)*f*r(1) + f*l^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);
            elseif id==2 % sel
                mu = mu + ((p-1)*f*l*r(1) + f*m*r(2) + (1-p*f)*l*r(3))*KEE;
                sig2 = sig2 + L*((p-1)*f*l^2*r(1) + f*m^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);
            elseif id==3 % 0
                mu = mu + ((p-1)*f*r(1) + f*r(2) + (1-p*f)*r(3))*KEE;
                sig2 = sig2 + L*((p-1)*f*r(1) + f*r(2) + (1-p*f)*r(3))*KEE*J(1,1);
            else
                rE = f*(r(2) + (p-1)*r(1)) + (1-p*f)*r(3);        
                mu = mu + tm(2)*eps(2,1)*N(1)*rE*J(2,1);
                sig2 = sig2 + L*tm(2)*eps(2,1)*N(1)*rE*J(2,1).^2;
            end

            % compute potential landscape 
            if id<4
                CModPar = num2cell(ModPar(1,:));
            else
                CModPar = num2cell(ModPar(2,:));
            end
            [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
            lb = min(EL,Vr) - 20.0;

            S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(id)/gL).*x + sf^2*exp((x-VT)./sf));
            q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
            dr(id) = (1000.0*sig2/(2*tau*q) - r(id))/(tau/1000.0); 
        end

        a = r2-r;
        fr(i) = dot(a,dr)/norm(a); % component along the FP-vector

        if i>1
            F = @(z) interp1(rates(1:i),fr(1:i),z,'spline');
            E(i) = -integral(@(x) F(x),rates(1),rates(i));
        else
            E(1) = 0;
        end
    end
    
elseif flg==1
    
    for i=1:length(rates)
        
        % compute fixed point for (N-1)-system
        r = zeros(1,4);
        qr = (rates(i)-FP(1,2))/(FP(2,2)-FP(1,2));   
        r0 = FP(1,:) + qr*(FP(2,:)-FP(1,:));
        
        r_opt = CompRate(ModPar,ConPar,LrnPar,r0([1 3 4]),rx,rates(i));
        r([1 3 4]) = r_opt;
        r(2) = rates(i);
        
        % influence of external pools
        mu = eps(1,NTypes)*N(NTypes)*J(1,NTypes)*rx(1)*tm(1);
        sig2 = J(1,NTypes)*mu*L;    

        % influence of inh pool
        mu = mu + tm(1)*eps(1,2)*N(2)*r(4)*J(1,2);
        sig2 = sig2 + L*tm(1)*eps(1,2)*N(2)*r(4)*J(1,2).^2;

        % influence of exc pools
        mu = mu + ((p-1)*f*l*r(1) + f*m*r(2) + (1-p*f)*l*r(3))*KEE;
        sig2 = sig2 + L*((p-1)*f*l^2*r(1) + f*m^2*r(2) + (1-p*f)*l^2*r(3))*KEE*J(1,1);

        % compute energy landscape along 1 dimension
        CModPar = num2cell(ModPar(1,:));
        [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
        lb = min(EL,Vr) - 20.0;
        
        S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*rates(i)/gL).*x + sf^2*exp((x-VT)./sf));
        q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
        fr(i) = 1000.0*sig2/(2*tau*q); 

        if i>1
            F = @(z) interp1(rates(1:i),fr(1:i),z,'spline');
            E(i) = -integral(@(x) (F(x)-x),rates(1),rates(i))/tm(1);
        else
            E(1) = 0;
        end

    end
    
end

end


%%%%%%%%%% helper functions %%%%%%%%%%

function rates = CompRate(ModPar,ConPar,LrnPar,r0,rx,r_fix)

options=optimset('TolFun',1e-5,'TolX',1e-5,'Display','off');
MSE = @(x) sum((CompTransferFun(x,rx,ModPar,ConPar,LrnPar,r_fix)-x).^2);
[rates,err] = fminsearch(MSE,r0,options);
if err>3*1e-10
    warning('Optimization error > 3*e-10')
end

end


function rates = CompTransferFun(r,rx,ModPar,ConPar,LrnPar,r_fix)

tm = ModPar(:,1)/1000.0;

f = LrnPar.f;
p = LrnPar.p;
m = LrnPar.m;
l = 1 - f*(m-1)/(1-f);
SDJ = LrnPar.SDJ;  
L = (1 + SDJ^2);

eps = ConPar.R;
N = ConPar.NN;
J = ConPar.J; 
NTypes = length(N);
rates = zeros(1,NTypes); 
KEE = tm(1)*eps(1,1)*N(1)*J(1,1);

% excitatory pools after learning (except for selective pool)
for i=1:2 % exc pools (+, 0)
    
    % influence of external pool
    mu = eps(1,NTypes)*N(NTypes)*J(1,NTypes)*rx(1)*tm(1);
    sig2 = J(1,NTypes)*mu*L;
    
    % influence of remaining populations
    for j=2:NTypes-1
        mu = mu + tm(1)*eps(1,j)*N(j)*r(1+j)*J(1,j);
        sig2 = sig2 + L*tm(1)*eps(1,j)*N(j)*r(1+j)*J(1,j).^2;
    end
    
    % influence of exc pools
    if i==1 % +
        mu = mu + (((p-2)*l + m)*f*r(1) + f*l*r_fix + (1-p*f)*l*r(2))*KEE;
        sig2 = sig2 + L*(((p-2)*l^2 + m^2)*f*r(1) + f*l^2*r_fix + (1-p*f)*l^2*r(2))*KEE*J(1,1);
    elseif i==2 % 0
        mu = mu + ((p-1)*f*r(1) + f*r_fix + (1-p*f)*r(2))*KEE;
        sig2 = sig2 + L*((p-1)*f*r(1) + f*r_fix + (1-p*f)*r(2))*KEE*J(1,1);
    end
     
    % compute transfer function  
    CModPar = num2cell(ModPar(1,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;

    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(i)/gL).*x + sf^2*exp((x-VT)./sf));
    q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
    rates(i) = 1000.0*sig2/(2*tau*q);

end
rE = f*(r_fix + (p-1)*r(1)) + (1-p*f)*r(2);

% remaining populations
for i=2:NTypes-1
    
    % influence of external pool
    mu = eps(i,NTypes)*N(NTypes)*J(i,NTypes)*rx(i)*tm(i);
    sig2 = J(i,NTypes)*mu*L;
    
    % influence of excitatory pools
    mu = mu + tm(i)*eps(i,1)*N(1)*rE*J(i,1);
    sig2 = sig2 + L*tm(i)*eps(i,1)*N(1)*rE*J(i,1).^2;
    
    % influence of remaining population
    for j=2:NTypes-1
        mu = mu + tm(i)*eps(i,j)*N(j)*r(1+j)*J(i,j);
        sig2 = sig2 + L*tm(i)*eps(i,j)*N(j)*r(1+j)*J(i,j).^2;
    end 
    
    % compute transfer function  
    CModPar = num2cell(ModPar(i,:));
    [tau,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
    lb = min(EL,Vr) - 20.0;

    S = @(x) (2.0/sig2)*(-(x.^2)/2 + (EL + mu - b*(tw/1000)*r(1+i)/gL).*x + sf^2*exp((x-VT)./sf));
    q = integral2(@(V,u) exp(S(V)-S(u)),lb,Vup,@(V) max(V,Vr),Vup);
    rates(1+i) = 1000.0*sig2/(2*tau*q);
end

end


