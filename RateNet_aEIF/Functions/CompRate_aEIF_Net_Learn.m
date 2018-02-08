function rates = CompRate_aEIF_Net_Learn(ModPar,ConPar,LrnPar,r0,rx,flag)
% This function computes the firing rate of an aEIF network after learning 
% (that is, several excitatory pools have been generated and respond to 
% specific stimuli) 
%
% ModPar = matrix of neuron model parameters (rows = neuron populations,
%          columns = parameters)
% ConPar = structure containing connection probabilities (ConPar.p), number
%          of neurons per population (ConPar.N) and 
%          connection strengths (ConPar.J)
% LrnPar = structure containing all parameters related to learning
% r0 = vector of initial values/first guesses
% rx = vector of external rates
% flag = factor variable (0: minimize MSE, 1: euler method)
% rates = vector of firing rates XXXX

tm = ModPar(:,1);
tm = [tm(1);tm(1);tm];
dt = min(tm)/10;
N = size(ModPar,1);

if flag==0
    options=optimset('TolFun',1e-5,'TolX',1e-5,'Display','off');
    MSE = @(x) sum((CompTransferFun_aEIF_Net_Learn(x,rx,ModPar,ConPar,LrnPar)-x).^2);
    [rates,err] = fminsearch(MSE,r0,options);
    if err>N*1e-10
        warning('Optimization error > Ne-10')
    end
elseif flag==1
    err=1.0;
    rates = r0;
    while(err>N*1e-10)
        f = CompTransferFun_aEIF_Net_Learn(rates,rx,ModPar,ConPar,LrnPar);
        rn = (1-dt./tm').*rates + dt*f./tm';
        err = sum((rn-rates).^2);
        rates = rn;
    end
else
    warning('Flag needs to be either 0 or 1.')
    return
end

end

