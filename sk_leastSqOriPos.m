function diffsq = sk_leastSqOriPos(eul_cntr,data2d,camParaCalib,param,model)

ncams = size(camParaCalib,2);
cg = eul_cntr(4:6);
d = param.diamet;
hl = param.armlength;
% imgd=zeros(ncams,1024,1280,'single');
% imgm=zeros(ncams,1024,1280,'single');

%make the model
eul = eul_cntr(1:3);
ori  = sk_ori(eul);

oarms = zeros(length(model.arms),3);
ocntr = zeros(length(model.cntr),3);

for i=1:length(model.arms)
    oarms(i,:) = (ori*model.arms(i,:)')';
    ocntr(i,:) = (ori*model.cntr(i,:)')';
end
oarms = bsxfun(@plus, oarms*hl,cg);
ocntr = bsxfun(@plus, ocntr*hl,cg);

diffsq = 0;

for icam=1:ncams
%     xpix = camParaCalib(icam).Npixw;
%     ypix = camParaCalib(icam).Npixh;
    % Determine the end points of the objects projected onto each camera
    oarms_proj = single((calibProj_Tsai(camParaCalib(icam),oarms))); % FLOAT
    ocntr_proj = single((calibProj_Tsai(camParaCalib(icam),ocntr)));

    % First we organize...
    % ... the experimental data
    
    pixd = single(data2d(icam).ind(:,1:2));
    vd = ones(size(data2d(icam).ind(:,3)));
    
    % We use the average intensity on each image to set the maximum value
    % in our model.
    I0 = single(mean(vd)); % FLOAT
%     I0 = single(1);
    % Compute the intensity of the model at the bright points in the data
    vm = sk_gaussian_intensity_multiple_rods(pixd,oarms_proj,ocntr_proj,I0,d);
    data2d(icam).indmod(:,3) = vm;
    
    plots = 0;
%     pos = zeros(size(vm));
%     pos(:) = sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
%     imgd(icam,pos) = data2d(icam).ind(:,3);
%     imgm(icam,pos)=vm;
    
    if plots
        center = (calibProj_Tsai(camParaCalib(icam),cg));
        figure(79);
        subplot(2,2,icam),imagesc(squeeze(imgd(icam,:,:)));
        colormap(gray);
        hold on;
%         plot(pixd(:,1),pixd(:,2),'wo','MarkerSize',1);
%         scatter(pixd(:,1),pixd(:,2),1,vd);
        plot([ocntr_proj(1,1);oarms_proj(1,1);ocntr_proj(2,1);oarms_proj(2,1);ocntr_proj(3,1);oarms_proj(3,1)],[ocntr_proj(1,2);oarms_proj(1,2);ocntr_proj(2,2);oarms_proj(2,2);ocntr_proj(3,2);oarms_proj(3,2)],'-w','LineWidth',2);
%         plot([x_L(:,1);x_R(end,1)],[x_L(:,2);x_R(end,2)],'-wo','LineWidth',2,'MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',7);
        axis([center(1)-128,center(1)+128,center(2)-160,center(2)+160]);
        set(gca,'FontSize',16);
        xlabel('x [pixel]');
        ylabel('y [pixel]');
        hold off;
    end
    
    % Now we compute the sum of the difference between the intesities
    % squared - at only the points just computed.
    
    diffsq = diffsq + sum((vm-vd).^2)*sqrt(size(data2d(icam).ind,1));
    
end

% pause(.1)

end
