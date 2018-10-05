%% PBPL PERiod AVErage 1D FEL simulation code %%
%%% Input deck intended to be compatible with WafFEL, 1D period average, and GENESIS %%%
%% P. Musumeci oscillator version %%
function [meanEff,Eff,blist] = perave_opti(R56buncher, phaseshift, cavitydetuning, transmission,npassmax)
% clear all
close all
GIT_dir;
figdir=[datadir,'perave_opti\'];
if exist(figdir,'dir')==0
    mkdir(figdir);
end


%% physical constants
global c mu0 e0 me eps0 IA Z0

c = 2.99792458.*10^8;                                            % speed of light
e0 = 1.60217657e-19;                                             % electron charge
me = 9.10938291e-31;                                             % electron mass
eps0 = 8.85418782e-12;                                           % eps_0
mu0 = 1.256637e-6;                                               % mu_0
IA = 17045;                                                      % Alfven current
Z0 = 376.73;                                                     % Impedance of free space         
power=[];
%% Load the User Determined initial conditions
clear power radfield thetap gammap bunch
param.sigma_t = 3e-13/2;
param.use3Dcorrection  = 1;
param.beamdistribution = 2;       % Using GENESIS flag: 2-uniform 1-gaussian
param.laserdistribution = 1;         % Using GENESIS flag: 2-uniform 1-gaussian
recirculate = 0;
t1 = tic;
Perave_User_Input_osc_pb;

%% Compute the undulator field
compute_undulator_field_v5h
firstpass=1;
loadtapering=0;
% perave_core_v6;
% perave_R56;
Perave_User_Input_osc;
% param.prebunching=1; %Turnoff prebuncher

compute_undulator_field_v5h

% param.nslices=size(oldfield,2);
firstpass=1;


%% Calculate 1-D FEL parameters
rho1D = 1/param.gamma0*(1/8*param.I/IA*param.K.^2/param.sigmax^2/param.ku^2)^(1/3);
Lgain = param.lambdau/(4*sqrt(3)*pi*rho1D);
Lsat =   param.lambdau/rho1D;
Psat = 1.6*rho1D*param.Ee*param.I;
if param.tapering
    [psi1, psi2, bucket_height, capture_fraction, bucket_area, bunching_factor] = bucket_parameters(param.psir);
    a1 = 2*param.lambda0/param.lambdau*e0*param.E0/me/c^2*sin(param.psir);
    a2 = ((2*param.lambda0/param.lambdau)^1.5)*Z0*e0*param.I*sin(param.psir)^2*capture_fraction*bunching_factor/2/param.A_e/me/c^2;
    pmax_prediction=P0+param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*param.Ee*param.I*capture_fraction;
    etamax = param.K*(a1*lwig+a2*lwig^2/2)/(1+param.K^2)*capture_fraction;
    bunchlength_rms = param.sigma_t;
    peakcurrent = param.I;
end
calculate_3Dcorrection; 



%% Run the main integration routine
% cavitydetuning = -30;    % In units of zsep
% transmission = 0.66;      % Power transmission through one cavity pass 
                                      % losses = 1 - transmission                                      
sigma_omega = 0.003*param.nslices*param.zsep;     % Filter fractional bandwidth. 


tapering_strength = 2;   % 0 max of slices at time 0 
                                      % 1 max of slices
                                      % 2 avg of slices

    

   updatetapering=0;
%     loadtapering=0;
%    load('D:\Matlab_data\TESSO_Kzload.mat');
 %% Oscillator loop
%  Kz_load=2.8202*ones(1,8);
%  loadtapering=0;
% npasses=1;
% perave_core_v6;
% perave_postprocessor_v6;
% 
% npassmax=60;
rad_vs_beam=zeros(param.nslices,npassmax);
getnewR56=0;
for npasses = 1:npassmax
%     clear power radfield thetap gammap bunch 
clear power radfield bunch
    t0 = tic;
    
Perave_User_Input_osc_pb;
if ~firstpass
param.nslices=size(thetap,2);
end


%% Compute the undulator field
compute_undulator_field_v5h
loadtapering=0;
% clear power radfielsd bunch

perave_core_v6;
perave_R56;
figure(7)
    title(titlestr);
    subplot(1,2,1)
        hold on
      plot(abs(fftshift(fft(oldfield))),'r');

  
        legend('oldfield')

    subplot(1,2,2)

    plot(abs(oldfield).^2/max(abs(oldfield).^2),'r')
            hold on

    plot(profile_b,'b')
    hold off
blist(npasses)=newbfactor;
Perave_User_Input_osc;
% param.prebunching=1; %Turnoff prebuncher
if ~firstpass
param.nslices=size(thetap,2);
end

compute_undulator_field_v5h

    perave_core_v6;
    disp(['Simulation time = ',num2str(toc(t0)./60),' min'])
    perave_postprocessor_v6   
    rad_vs_und(:,npasses) = sum(power,2)*param.lambda0*param.zsep/c;
    
    rad_vs_beam(:,npasses) = power(end,:);
    Eff(npasses) = Efficiency;
    PL(npasses) = pulselength;
    oldfield=zeros(1,param.nslices);
    if param.itdp==1
    if  cavitydetuning>0
    oldfield(1,cavitydetuning+1:cavitydetuning+size(radfield,2)) = radfield(end,:)*sqrt(transmission);
    else
    oldfield(1,1:1+cavitydetuning+size(radfield,2)) = radfield(end,-cavitydetuning:end)*sqrt(transmission);    
    end
    else
        oldfield=radfield(end,:)*sqrt(transmission);
    end
    pause(0.5)
    
    
    perave_filter;

    %%
    figure(10)
    titlestr=sprintf('npass=%.f cavitydetuning=%.2f transmission=%.2f Q=%.2e',npasses,cavitydetuning,transmission,Q);
    title(titlestr);
    subplot(1,2,1)
        hold on
      plot(abs(fftshift(fft(oldfield))),'r');
    plot(abs(fftshift(fft(oldfield)).*filter3),'g');

  
        legend('oldfield','filterfield')

    subplot(1,2,2)
    
    filterfield = ifft(ifftshift(fftshift(fft(oldfield) ).*filter3));
    plot(power(end,:),'k')
    hold on
    plot(abs(filterfield).^2/377*param.A_e,'g')
    plot(abs(oldfield).^2/377*param.A_e,'r')
    plot(profile_b*max(power(end,:))*0.5,'b')
    hold off
    pause(0.5)
    legend('power','filterfield','oldfield','profile_b')
    if param.itdp==1
    oldfield = filterfield;
    end
    firstpass = 0;                                  % Start recirculation
%         saveas(gcf,[figdir,'field_',num2str(npasses),'.png'])
%         figure(2)
% %                 saveas(gcf,[figdir,'outfig_',num2str(npasses),'.png'])
% 
%         figure(3)
% %                         saveas(gcf,[figdir,'contour_',num2str(npasses),'.png'])
% 
%         figure(4)
% %                                 saveas(gcf,[figdir,'spec_',num2str(npasses),'.png'])
                                
%                                 if npasses>1 & abs((mean(rad_vs_beam(:,npasses))-mean(rad_vs_beam(:,npasses-1)))/mean(rad_vs_beam(:,npasses)))<.01
%                                     updatetapering=1;
%                                 else
%                                     updatetapering=-1;
%                                 end
%                                 
%                                 if Eff>0.25
%                                     updatetapering=-1;
%                                 end

% if npasses>1 & npasses<50
%     updatetapering=-1;
% elseif npasses==50
%     updatetapering=1;
% else 
%     updatetapering=-1;
% end
    
%                                     
Kz_save(:,npasses)=Kz;
% if npasses>=20
% getnewR56=1;
% end

end
perave_opti_str=sprintf('Effmean=%.2e EffEnd=%.2e bfend=%.2f bfmean=%.2f r56=%.2e ps=%.2f cd=%.f trans=%.2f ',mean(Eff),Eff(end),newbfactor,mean(blist),R56buncher,phaseshift,cavitydetuning,transmission);

%% Post-process stuff
figure(100)
plot(max(rad_vs_und),'b')
title('max rad vs und')
        saveas(gcf,[figdir,perave_opti_str,'final_radvsund.png'])

figure(101)
plot([1:1:param.Nsnap]*param.stepsize,rad_vs_und(:,end),'r')
hold on
plot([1:1:param.Nsnap]*param.stepsize, meanenergy*charge*511000)
xlim([0,param.Nsnap*param.stepsize])
title('Radiation energy along undulator')
        saveas(gcf,[figdir,perave_opti_str,'final_radenergy.png'])
hold off
figure(102)
plot(PL)
title('pulselength')
        saveas(gcf,[figdir,perave_opti_str,'final_pulselength.png'])

figure(103)
plot(Eff)
title('Eff')

        saveas(gcf,[figdir,perave_opti_str,'final_eff.png'])
if param.itdp
figure(300)
contourf([1:size(rad_vs_beam,1)]*param.zsep*param.lambda0/c,[1:npassmax],rad_vs_beam')
title('rad vs beam')
        saveas(gcf,[figdir,perave_opti_str,'final_beam.png'])

colorscheme=cool(size(rad_vs_und,2));
hold off

end

figure(104)
plot(blist)
hold off
title('bunch factor in each run')
        saveas(gcf,[figdir,perave_opti_str,'bfactor.png'])
meanEff=-mean(Eff);
meanEff=-blist(end);
end