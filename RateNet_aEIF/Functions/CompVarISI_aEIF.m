function VarISI = CompVarISI_aEIF(r,ModPar,InpPar)
% This function computes the variance of the ISIs for a aEIF neuron in the
% presence of Gaussian white noise
%
% r = mean firing rate (1/r = mean ISI)
% ModPar = vector of neuron model parameters
% InpPar = vector of input mean and standard deviation
% VarISI = variance of ISIs

% rename properly
CModPar = num2cell(ModPar);
CInpPar = num2cell(InpPar);
[tm,gL,EL,sf,VT,Vr,Vup,b,tw] = CModPar{:};
[mu,sig] = CInpPar{:};

% compute averaged adaptation & E & lb
avgw = b*r*tw/1000.0;
E = EL + mu -avgw/gL;
lb = min(EL,Vr) - 20.0;

% compute ISI-variance
ATol=1e-4;
RTol=1e-2;

VarISI = integral2(@InnerFunc,Vr,Vup,lb,@(V) V,'Method','tiled','AbsTol',ATol,'RelTol',RTol);

    function out = InnerFunc(V,u)  

        F=@(r) (2*E.*r - r.^2 + 2*sf^2.*exp((r-VT)./sf))/sig^2;
        I=arrayfun(@(u) integral(@(z) exp(F(z)-F(u)), lb, u, 'ArrayValued', true, 'AbsTol',ATol,'RelTol',RTol),u);

        out = (8*(tm/1000.0)^2*exp(F(u)-F(V)).*I.^2)/sig^4;
    end

 end

