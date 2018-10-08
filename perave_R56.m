thetap_R56=[];
drift_buncher=28*.032;
%%change R56?
% if getnewR56==1
%                              Klz = mean(abs(radfield(end,1:param.nslices)),2);
%                              Klz_mod=Klz*pi*param.lambda0/2*param.chi2;
%                              Klz_mod=Klz/2/pi*param.lambda0*param.chi2;
%     param.R56buncher=(1/pi)*param.gamma0^2*param.lambda0^2/(2*pi)/(param.K*(lwig)*Klz_mod)+drift_buncher/(2*param.gamma0^2);
% R56buncher=0;
% else
%     R56buncher=param.R56buncher;
% end

for islice = 1:size(gammap,2)
    gamma_avg=mean(gammap(end,islice,:));
for ip = 1:param.Np
    
    thetap_R56(islice,ip)=thetap(end,islice,ip)+  2*pi*param.R56buncher/param.lambda0*(gammap(end,islice,ip)-gamma_avg)/gamma_avg+param.phaseshift;        
end
    bcomplex_old(islice) = (sum(exp(1i.*thetap(end,islice,:))/param.Np));
    bcomplex_new(islice) = (sum(exp(1i.*thetap_R56(islice,:))/param.Np));

end

gammapend=gammap(end,:,:);
thetapend=thetap(end,:,:);
gammapend=reshape(gammapend,[size(gammap,2),param.Np]);
thetapend=reshape(thetapend,[size(gammap,2),param.Np]);
thetap_R56=reshape(thetap_R56,[size(gammap,2),param.Np]);
figure(88)
plot(thetapend(round(islice/2),:),gammapend(round(islice/2),:))
hold on
plot(thetap_R56(round(islice/2),:),gammapend(round(islice/2),:))
thetap(end,:,:)=thetap_R56;
oldbfactor=mean(abs(bcomplex_old));
newbfactor=mean(abs(bcomplex_new));

blist(npasses)=newbfactor;

titlestr=sprintf('R56buncher=%.2e oldbunching=%.2f newbunching=%.2f',param.R56buncher,oldbfactor,newbfactor);
title(titlestr);
hold off


%%%slippage due to prebuncher
if param.itdp==1
radfield_R56=[];
R56slippage = round(param.R56buncher/param.lambda0)*2/param.zsep;
radfield_R56(1,1:R56slippage-1)=0;    
radfield_R56(1,R56slippage:size(radfield,2)) = radfield(end,1:(size(radfield,2)-R56slippage+1));

figure(9)
plot(abs(radfield(end,:)));
hold on
plot(abs(radfield_R56));
hold off

    oldfield=radfield_R56;
      
else
    oldfield=radfield(end,:);

end

param.prebunched=1;
    
%     thetap(:,1:R56slippage,:)=[];
%     gammap(:,1:R56slippage,:)=[];
%     thetap_R56(1:R56slippage,:)=[];
%     gammapend(1:R56slippage,:)=[];
%     bunch(:,1:R56slippage)=[];
%     profile_b(1:R56slippage)=[];
%     profile_l(1:R56slippage)=[];
    