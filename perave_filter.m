
%% Filter definition (Filter2 is a complex transfer function. Cavity detuning needs to be adjusted to 12)
    jfreq = 1:param.nslices;
    filter = exp(-(jfreq-param.nslices/2).^2/2/sigma_omega^2);
    filter2=[];filter3=[];
    for jfreq = 1:param.nslices
    y = (jfreq-param.nslices/2)/sigma_omega;
    if(y>=1)
        filter2(jfreq) = y-sqrt(y.^2-1); %ryan lindberg (KJ Kim) bragg mirror
    elseif(y<=-1)
        filter2(jfreq) = (y+sqrt(y.^2-1));
    else
        filter2(jfreq) = y+1i*sqrt(1-y.^2);
    end
        omega_m=param.nslices/2;
        Q = 0.3;
        filter3(jfreq) = 1i*jfreq/Q / (omega_m^2-jfreq^2+1i*jfreq/Q);   %dispersion
    end

    
    filterdelay = round(param.nslices/2/pi/sigma_omega);
    figure(200)
    plot(filter)
    hold on
    plot(abs(filter2))
    plot(abs(filter3),'k')
    hold off
    legend('filter','filter2','filter3')
    GIT_dir
    figure(201)
    plot(angle(filter2))
    hold on
    plot(angle(filter3),'k')
    hold off