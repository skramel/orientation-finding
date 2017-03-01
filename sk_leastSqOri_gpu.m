function diffsq = sk_leastSqOri_gpu(eul,cg,data2d,camParaCalib,model,ICS_treshold)

ncams = size(camParaCalib,2);

img = zeros(ncams,1024,1280,'single');

%make the model
ori  = sk_ori(eul);

arms_ori = zeros(size(model.arms));
cntr_ori = zeros(size(model.cntr));

d = model.rad;
    
for i=1:size(model.arms,1)
    arms_ori(i,:) = (ori*model.arms(i,:)')';
%     cntr_ori(i,:) = (ori*model.cntr(i,:)')';
end
arms = bsxfun(@plus, arms_ori,cg);
cntr = bsxfun(@plus, cntr_ori,cg);
mid = (cntr+arms)./2;

diffsq = 0;

for icam=1:ncams

    % Determine the end points of the objects projected onto each camera
    arms_proj = single((calibProj_Tsai(camParaCalib(icam),arms)));
    cntr_proj = single((calibProj_Tsai(camParaCalib(icam),cntr)));
    
    % Determine the radius of the objects projected onto each camera
    rad = cross(arms_ori,repmat(camParaCalib(icam).Tinv',3,1),2);
    rad = bsxfun(@rdivide,rad,sqrt(sum(rad.^2,2))).*d;
    rad = bsxfun(@plus, rad, mid);
    rad_proj = single((calibProj_Tsai(camParaCalib(icam),rad)));
    
    % First we organize the experimental data
    pixd = single( data2d(icam).ind(:,1:2) );
    vd = single( data2d(icam).ind(:,3) );
%     vd = ones(size(data2d(icam).ind(:,3)));
    
    % We use the average intensity on each image to set the maximum value
    % in our model.
%     I0 = single( ICS_treshold(icam)-mean(data2d(icam).ind(:,3)) );
    I0 = single( ICS_treshold(icam)-90 );
%     I0 = single(1);
    % Compute the intensity of the model at the bright points in the data
    vm = sk_fermi_intensity_multiple_rods_gpu(pixd, arms_proj, rad_proj, cntr_proj, I0, ICS_treshold(icam));
%     data2d(icam).indmod(:,3) = vm;

    % Now we compute the sum of the difference between the intesities
    % squared - at only the points just computed, weighted by the number of
    % bright pixel of the data on each camera
    diffsq = diffsq + sum((vd-vm).^2);
    
%     % plotting
%     xpix = camParaCalib(icam).Npixw;
%     ypix = camParaCalib(icam).Npixh;
%     pos = sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
% %     img(icam,pos) = data2d(icam).ind(:,3);
%     img(icam,pos) = vm;
%     
%     center = (calibProj_Tsai(camParaCalib(icam),cg));
%     figure(79);
%     subplot(2,2,icam),imagesc(squeeze(img(icam,:,:)));
%     colormap(gray);
%     hold on;
% 
%     plot([cntr_proj(1,1);arms_proj(1,1)],[cntr_proj(1,2);arms_proj(1,2)],'-w','LineWidth',2);
%     plot([cntr_proj(2,1);arms_proj(2,1)],[cntr_proj(2,2);arms_proj(2,2)],'-w','LineWidth',2);
%     plot([cntr_proj(3,1);arms_proj(3,1)],[cntr_proj(3,2);arms_proj(3,2)],'-w','LineWidth',2);
% 
%     axis([center(1)-256,center(1)+256,center(2)-320,center(2)+320]);
%     set(gca,'FontSize',16);
%     xlabel('x [pixel]');
%     ylabel('y [pixel]');
%     hold off;
%     % end plotting
    
end

diffsq = gather(diffsq);
% pause(.1) % don't forget to comment out if not plotting!

end
