thetap_new=[];
for islice = 1:size(gammap,2)
    gamma_avg=mean(gammap(end,islice,:));
for ip = 1:param.Np
        thetap_new(islice,ip)=thetap(end,islice,ip)+  2*pi*param.R56buncher/param.lambda0*(gammap(end,islice,ip)-gamma_avg)/gamma_avg+param.phaseshift;        
end
bcomplex_old(islice) = (sum(exp(1i.*thetap(end,islice,:))/param.Np));
bcomplex_new(islice) = (sum(exp(1i.*thetap_new(islice,:))/param.Np));

end

gammapend=gammap(end,:,:);
thetapend=thetap(end,:,:);
gammapend=reshape(gammapend,[size(gammap,2),param.Np]);
thetapend=reshape(thetapend,[size(gammap,2),param.Np]);
thetap_new=reshape(thetap_new,[size(gammap,2),param.Np]);
figure(8)
plot(thetapend(round(islice/2),:),gammapend(round(islice/2),:))
hold on
plot(thetap_new(round(islice/2),:),gammapend(round(islice/2),:))
thetap(end,:,:)=thetap_new;
oldbfactor=mean(abs(bcomplex_old));
newbfactor=mean(abs(bcomplex_new));
titlestr=sprintf('oldbunching=%.2f newbunching=%.2f',oldbfactor,newbfactor);
title(titlestr);


%%%slippage in prebuncher
radfield_new=[];
R56slippage = round(param.R56buncher/param.lambda0)*2/param.zsep;
    radfield_new(1,1:R56slippage-1)=0;
%     if firstpass
    radfield_new(1,R56slippage:size(radfield,2)) = radfield(end,1:(size(radfield,2)-R56slippage+1));
%     else
%             radfield_new(1,R56slippage:size(oldfield,2)) = oldfield(end,1:(size(oldfield,2)-R56slippage+1));
%     end
    figure(9)
    plot(abs(radfield(end,:)));
    hold on
    plot(abs(radfield_new));
    hold off
%     radfield_new(:,1:R56slippage)=[];

    oldfield=radfield_new;
    firstpass=0;
   
    param.nslices=size(oldfield,2);
%     
%     thetap(:,1:R56slippage,:)=[];
%     gammap(:,1:R56slippage,:)=[];
%     thetap_new(1:R56slippage,:)=[];
%     gammapend(1:R56slippage,:)=[];
%     bunch(:,1:R56slippage)=[];
%     profile_b(1:R56slippage)=[];
%     profile_l(1:R56slippage)=[];
%     