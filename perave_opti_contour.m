transmission_list=.01:.02:.99;
R56_list=1:.5:30;

for j1=1:1:length(R56_list)
    for j2=1:1:length(transmission_list)
[Effend,Efflist,blist]=perave_opti(R56_list(j1)*1e-6,1.8,-20,transmission_list(j2),30);
Eff_totlist(j1,j2,:)=Efflist;
blist_totlist(j1,j2,:)=blist;
Eff_end_totlist(j1,j2)=Efflist(end);
blist_end_totlist(j1,j2)=blist(end);
Eff_mean_totlist(j1,j2)=mean(Efflist);
blist_mean_totlist(j1,j2)=mean(blist);
    end
end

figure
contourf(R56_list,transmission_list,Eff_end_totlist);

figure
contourf(R56_list,transmission_list,blist_end_totlist);