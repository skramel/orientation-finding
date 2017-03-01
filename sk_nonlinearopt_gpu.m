function [eul_cntr,fval,out,xflag] = sk_nonlinearopt_gpu(eul_cntr,data2d,camParaCalib,model,ICS_treshold,check)

    eul = eul_cntr(1:3);
    cntr = eul_cntr(4:6);

    TolX = 1e-4;
    TolFun = 1e-4;
    nonlin = tic;
    
    options = optimset('TolX',TolX,'TolFun',TolFun,'Display','final','LargeScale','off'); 
    % [eul,fval1,~,~] = fminsearch(@(e) sk_leastPixOri_gpu(e,data2d,camParaCalib,model,cntr,ICS_treshold), eul, options);
    [eul,fval1,~,~] = fminsearch(@(e) sk_leastSqOri_gpu(e,cntr,data2d,camParaCalib,model,ICS_treshold), eul, options);
%     if check
%         eul2 = rand(1,3)*2*pi;
%         % [eul2,fval2,~,~] = fminsearch(@(e) sk_leastPixOri_gpu(e,data2d,camParaCalib,model,cntr,ICS_treshold), eul2,options);
%         [eul2,fval2,~,~] = fminsearch(@(e) sk_leastSqOri_gpu(e,cntr,data2d,camParaCalib,model,ICS_treshold), eul2, options);         
%         if fval2<fval1
%             eul_cntr = [eul2 cntr];
%         else
%             eul_cntr = [eul cntr];
%         end
%     else
%         eul_cntr = [eul cntr];
%     end
        
    display(sprintf('\tElapsed time for first fminsearch: %f sec', toc(nonlin)));
    
    TolX = 1e-6;
    TolFun = 1e-6;    
    nonlin2 = tic;
    
    options = optimset('TolX',TolX,'TolFun',TolFun,'Display','final','LargeScale','off'); 
    [eul_cntr,fval,xflag,out] = fminsearch(@(e) sk_leastPixOriPos_gpu(e,data2d,camParaCalib,model,ICS_treshold), eul_cntr,options);
    display(sprintf('\tElapsed time for second fminsearch: %f sec', toc(nonlin2)));
%     sk_plot_model(eul_cntr,data2d,camParaCalib,param,model);
    
    eul = eul_cntr(1:3);
%     Checking to make sure no Euler angles are outside of [-pi,pi]
    for i = 1:3
           if eul(i) > 2*pi
               eul(i) = mod(eul(i),2*pi);
           elseif eul(i) < -2*pi
               eul(i) = mod(eul(i),-2*pi);
           end
    end

end
