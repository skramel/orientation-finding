function int_out = sk_intensity_gpu(arm_co,armw,arml,ICS_treshold, I0)

    int_out = ICS_treshold-I0*(1./(1+exp((arm_co(:,2)-armw)./(0.1*armw)))).*(1./(1+exp((arm_co(:,1)-.5*arml)./(0.02*arml))));
