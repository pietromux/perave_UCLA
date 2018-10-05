function [newphasespace,newevalue]=push_FEL_particles_RK4(phasespace,evalue,param,kvalue)
   
gammar_sq=param.lambdau/(2*param.lambda0)*(1+kvalue^2);
sc = 1;
% Euler method for field (for some reason the most accurate...)
newevalue=evalue-param.stepsize*(param.chi1*kvalue*...
sc*mean(exp(-1i*phasespace(:,1))./phasespace(:,2)));

% Leapfrog method for the field


% RK-2 for the field
%k1e=-1*param.chi1*kvalue*mean(exp(-1j*phasespace(:,1))./phasespace(:,2));
%y1e=evalue+k1e*param.stepsize/2;
%k2e=-1*param.chi1*kvalue*mean(exp(-1j*(phasespace(:,1)+param.stepsize/2))./(phasespace(:,2)+param.stepsize/2));
%newevalue=evalue+k2e*param.stepsize;

% RK-4 for the particles

k1theta=param.stepsize*(param.ku*(1-(gammar_sq./phasespace(:,2).^2)));
k1gamma=param.stepsize*(param.chi2*(kvalue./phasespace(:,2)).*...
    real(evalue*sc*exp(1i*phasespace(:,1))));

k2theta=param.stepsize*(param.ku*(1-(gammar_sq./(phasespace(:,2)+0.5*k1gamma).^2)));
k2gamma=param.stepsize*(param.chi2*(kvalue./(phasespace(:,2)+0.5*k1gamma)).*...
    real(evalue*sc*exp(1i*(phasespace(:,1)+0.5*k1theta))));

k3theta=param.stepsize*(param.ku*(1-(gammar_sq./(phasespace(:,2)+0.5*k2gamma).^2)));
k3gamma=param.stepsize*(param.chi2*(kvalue./(phasespace(:,2)+0.5*k2gamma)).*...
    real(evalue*sc*exp(1i*(phasespace(:,1)+0.5*k2theta))));

k4theta=param.stepsize*(param.ku*(1-(gammar_sq./(phasespace(:,2)+k3gamma).^2)));
k4gamma=param.stepsize*(param.chi2*(kvalue./(phasespace(:,2)+k3gamma)).*...
    real(evalue*sc*exp(1i*(phasespace(:,1)+k3theta))));

newphasespace(:,1)=phasespace(:,1)+1/6*(k1theta+2*k2theta+2*k3theta+k4theta);
newphasespace(:,2)=phasespace(:,2)+1/6*(k1gamma+2*k2gamma+2*k3gamma+k4gamma);

% RK-4 for the field
% k1e=-1*param.stepsize*param.chi1*kvalue*mean(exp(-1j*phasespace(:,1))./phasespace(:,2));
% k2e=-1*param.stepsize*param.chi1*kvalue*mean(exp(-1j*(phasespace(:,1)+0.5*param.stepsize))./(phasespace(:,2)+0.5*param.stepsize));
% k3e=-1*param.stepsize*param.chi1*kvalue*mean(exp(-1j*(phasespace(:,1)+0.5*param.stepsize))./(phasespace(:,2)+0.5*param.stepsize));
% k4e=-1*param.stepsize*param.chi1*kvalue*mean(exp(-1j*(phasespace(:,1)+param.stepsize))./(phasespace(:,2)+param.stepsize));
% 
% newevalue=evalue+1/6*(k1e+2*k2e+2*k3e+k4e);
 