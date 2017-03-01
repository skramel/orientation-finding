function diffsq = sk_leastSqOriPos_triad(eul_cntr,data2d,camParaCalib,model)
ncams = size(camParaCalib,2);

%make the model
eul = eul_cntr(1:3);
cntr = eul_cntr(4:6);
ori  = sk_ori(eul);

d = model.rad;

triad_ori = (ori*model.ends')';
triad = triad_ori+repmat(cntr,size(triad_ori,1),1);
mid = repmat(cntr,size(triad_ori,1),1);

diffsq = 0;

for icam=1:ncams

    % Determine the end points and radius of the objects projected onto each camera
    rad = cross(triad_ori,repmat(camParaCalib(icam).Tinv',size(triad,1),1),2);
    rad = bsxfun(@rdivide,rad,sqrt(sum(rad.^2,2))).*d;
    rad = bsxfun(@plus, rad, mid);
    points_proj = single((sk_calibProj(camParaCalib(icam),[triad;mid;rad])));
    
    % First we organize the experimental data
    pixd = single( data2d(icam).ind(:,1:2) );
    vd = single( data2d(icam).ind(:,3) );
%     vd = ones(size(data2d(icam).ind(:,3)));
    
    I0 = single( mean(data2d(icam).ind(:,3)) );
%     I0 = single(1);
    % Compute the intensity of the model at the bright points in the data
    vm = sk_fermi_intensity_triad(pixd,points_proj,I0);
    data2d(icam).indmod(:,3) = vm;
        
    % Now we compute the sum of the difference between the intesities
    % squared - at only the points just computed, weighted by the number of
    % bright pixel of the data on each camera
    diffsq = diffsq + sum((vd-vm).^2);

end

end
