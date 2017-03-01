function diffPix = sk_leastPixOriPos_gpu(eul_cntr,data2d,camParaCalib,model,ICS_treshold)

ncams = size(camParaCalib,2);

%make the model
eul = eul_cntr(1:3);
cg = eul_cntr(4:6);

%make the model
ori  = sk_ori(eul);

ends = zeros(size(model.ends));
begs = zeros(size(model.begs));
narms = size(model.ends,3);

for i=1:narms
    ends(:,:,i) = bsxfun(@plus, (ori*model.ends(:,:,i)')', cg);
    begs(:,:,i) = bsxfun(@plus, (ori*model.begs(:,:,i)')', cg);
end

arms = [begs;ends];
    
diffPix = 0;

for icam=1:ncams
    
    xpix = camParaCalib(icam).Npixw;
    ypix = camParaCalib(icam).Npixh;
    % distance from origin along line of sight of camera
%     depth = abs(dot(camParaCalib(icam).Tinv./norm(camParaCalib(icam).Tinv),cg))./25; % divide by constant scale factor for gaussian blur
    
    % First we organize the experimental data
    pixd = single(data2d(icam).ind(:,1:2));
    pixd_ind = sub2ind([ypix,xpix],pixd(:,2),pixd(:,1));
    
    % Determine modelarms projected onto each camera
    for i=1:narms
        arms_proj = calibProj_Tsai(camParaCalib(icam),squeeze(arms(:,:,i)));
        DT = delaunayTriangulation(arms_proj);
        K = convexHull(DT);
        pixm = single(poly2mask(DT.Points(K,1),DT.Points(K,2),ypix,xpix));
%         pixm = imgaussfilt(pixm,depth);
        if i==1
            pixm_ind = find(pixm>0);
        else
            pixm_ind = cat(1,pixm_ind,find(pixm>0));
        end
    end
    
    img = zeros(1024,1280);
    img(pixd_ind)=single(data2d(icam).ind(:,3));
    img(pixm_ind)=100;
    subplot(2,2,icam);
    imagesc(img);colormap(gray);
    % Now we compute the sum of the difference between the pixel positions
    diffPix = diffPix + size(setdiff(pixd_ind,pixm_ind),1);
    
end

pause(.1)

end
