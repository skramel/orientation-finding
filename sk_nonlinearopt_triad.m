function [eul_cntr,fval,out,xflag] = sk_nonlinearopt_triad(eul_cntr,data2d,camParaCalib,model,check)

    eul = eul_cntr(1:3);
    cntr = eul_cntr(4:6);
    
    if check
        
        TolX = 1e-3;
        TolFun = 1e-3;
        MaxFunEvals = 250;
        nonlin = tic;

%         options = optimset('TolX',TolX,'TolFun',TolFun,'Display','final','LargeScale','off');
        options = optimset('TolX',TolX,'TolFun',TolFun,'MaxFunEvals',MaxFunEvals,'Display','final');
        [eul,fval1,~,~] = fminsearch(@(e) sk_leastSqOri_triad(e,cntr,data2d,camParaCalib,model), eul, options);
    
        eul2 = rand(1,3)*2*pi;
        [eul2,fval2,~,~] = fminsearch(@(e) sk_leastSqOri_triad(e,cntr,data2d,camParaCalib,model), eul2, options);
        if fval2<fval1
            eul_cntr = [eul2 cntr];
        else
            eul_cntr = [eul cntr];
        end
        
        display(sprintf('\tElapsed time for first fminsearch: %f sec', toc(nonlin)));
    else
        eul_cntr = [eul cntr];
    end

    TolX = 1e-3;
    TolFun = 1e-3;    
    nonlin2 = tic;
    
    options = optimset('TolX',TolX,'TolFun',TolFun,'Display','final'); 
    [eul_cntr,fval,xflag,out] = fminsearch(@(e) sk_leastSqOriPos7_triad(e,data2d,camParaCalib,model), eul_cntr,options);
    display(sprintf('\tElapsed time for second fminsearch: %f sec', toc(nonlin2)));
    
%     eul = eul_cntr(1:3);
%     cntr = eul_cntr(4:6);

%     Checking to make sure no Euler angles are outside of [-pi,pi]
%     for i = 1:3
%            if eul(i) > 2*pi
%                eul(i) = mod(eul(i),2*pi);
%            elseif eul(i) < -2*pi
%                eul(i) = mod(eul(i),-2*pi);
%            end
%     end

% % plotting
% 
%     ncams = size(camParaCalib,2);
%     img = zeros(ncams,1024,1280,'single');
% 
%     d = model.rad;
%     
%     ori  = sk_ori(eul);
%     
%     triad_ori = (ori*model.ends')';
%     triad = triad_ori+repmat(cntr,size(triad_ori,1),1);
%     mid = repmat(cntr,size(triad_ori,1),1);
% 
%     for icam=1:ncams
% 
%         xpix = camParaCalib(icam).Npixw;
%         ypix = camParaCalib(icam).Npixh;
% 
%         % Determine the end points and radius of the objects projected onto each camera
%         rad = cross(triad_ori,repmat(camParaCalib(icam).Tinv',size(triad,1),1),2);
%         rad = bsxfun(@rdivide,rad,sqrt(sum(rad.^2,2))).*d;
%         rad = bsxfun(@plus, rad, mid);
%         points_proj = single((calibProj_Tsai(camParaCalib(icam),[triad;mid;rad])));
%         
%         pixd = single( data2d(icam).ind(:,1:2) );
%         vd = single( data2d(icam).ind(:,3) );
% 
% %         I0 = single( mean(data2d(icam).ind(:,3)) );
% %         vm = sk_fermi_intensity_multiple_fiber3(pixd,points_proj,I0);
%         
%         pos = sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
%         img(icam,pos) = vd;
% 
%         figure(79);
%         subplot(2,2,icam),imagesc(squeeze(img(icam,:,:)));
%         colormap(gray);
%         hold on;
% %         points_proj(1:3,:) arms
% %         points_proj(4:6,:) cntr
% %         points_proj(7:9,:) rad
%         plot(points_proj(1:3,1),points_proj(1:3,2),'or','MarkerSize',7);
%         plot(points_proj(7:9,1),points_proj(7:9,2),'xr','MarkerSize',5);
% 
%         plot([points_proj(4,1);points_proj(1,1)],[points_proj(4,2);points_proj(1,2)],'-w','LineWidth',2);
%         plot([points_proj(5,1);points_proj(2,1)],[points_proj(5,2);points_proj(2,2)],'-w','LineWidth',2);
%         plot([points_proj(6,1);points_proj(3,1)],[points_proj(6,2);points_proj(3,2)],'-w','LineWidth',2);
% 
%         axis([points_proj(4,1)-96,points_proj(4,1)+96,points_proj(4,2)-120,points_proj(4,2)+120]);
% 
%         set(gcf, 'Position', get(0, 'Screensize'));
%         set(gca,'FontSize',16);
%         xlabel('x [pixel]');
%         ylabel('y [pixel]');
%         hold off;
%     end
%     pause(.1);
% 
% % end plotting

end
