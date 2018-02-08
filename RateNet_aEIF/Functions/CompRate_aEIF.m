function rate = CompRate_aEIF(ModPar,InpPar,r0,flag)
% This function computes the firing rate of an aEIF neuron in the presence
% of Gaussian white noise self-consistently
%
% ModPar = vector of neuron model parameters
% InpPar = vector of input mean and standard deviation
% r0 = initial value/first guess
% flag = factor variable (0: minimize MSE, 1: euler method)
% rate = firing rate of the neuron

tm = ModPar(1);
dt = tm/10;

if flag==0
    options=optimset('TolFun',1e-5,'TolX',1e-5,'Display','off');
    MSE = @(x) (CompTransferFun_aEIF(x,ModPar,InpPar)-x)^2;
    [rate,~] = fminsearch(MSE,r0,options);
elseif flag==1
    err=1.0;
    rate = r0;
    while(err>1e-10)
        f = CompTransferFun_aEIF(rate,ModPar,InpPar);
        rn = (1-dt/tm)*rate + dt*f/tm;
        err = (rn-rate).^2;
        rate = rn;
    end
else
    warning('Flag needs to be either 0 or 1.')
    return
end

end

