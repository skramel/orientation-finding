clear;close all;restoredefaultpath;

NAS = 'U';
addpath(genpath(['Z:\skramel\codes\'])); % codes
addpath([NAS ':\skramel\9-8-2016\results\']); % data

% define paths
filepath = [NAS ':\skramel\9-8-2016\results\'];
savepath = [NAS ':\skramel\9-8-2016\results\'];

filename = 'run1';

% Load in all the cpv files.
cpv1 = fopen([NAS ':\skramel\9-8-2016\9-8-2016-cam1\' filename '.cpv'],'r');
cpv2 = fopen([NAS ':\skramel\9-8-2016\9-8-2016-cam2\' filename '.cpv'],'r');
cpv3 = fopen([NAS ':\skramel\9-8-2016\9-8-2016-cam3\' filename '.cpv'],'r');
cpv4 = fopen([NAS ':\skramel\9-8-2016\9-8-2016-cam4\' filename '.cpv'],'r');

cpvs = [cpv1, cpv2, cpv3, cpv4];
ncams = length(cpvs);

% Calibration file
camParaCalib = load('Z:\skramel\codes\9-8-2016\calibration\camParaCalib-9-8-2016.mat');
camParaCalib = camParaCalib.camParaCalib;
% Define the number of pixels in each image. We assume that it is the same
% for each camera, which it currently is. (2012)
xpix = camParaCalib(1).Npixw;
ypix = camParaCalib(1).Npixh;

cpv_index_good = load([filepath filename '_cpv_index_good.mat']);
cpv_index_good = cpv_index_good.cpv_index_good;

% set parameters
ICS_treshold = [200 200 200 200]; % Image compression system threshold. This sets background brightness for model and data
cls_min = [1000 1000 1000 1000]; % min num of bright pixel for cluster
cls_max = [2500 2500 2500 2500]; % max num of bright pixel for cluster
cls_treshold = 1000; % min number of bright pixels required on each camera
cls_sep = 50;

fps = 450;

% create model
armlength = 9; % armlength of particle in mm
armdiam = .9; % arm diameter of particle in mm (if the diameter is kept a little larger and the "chemical potential" in the fermi function larger (softer edge) , the orientation finding codes work better)
% facets = 10;
% model = sk_triad4(armlength,armdiam,facets);
model = sk_triad3(armlength,armdiam);

data = read2_gdf([NAS ':\skramel\9-8-2016\results\' filename '\' filename '_clusters_final_dyn2_tracked_4cams.gdf']);
data(data(:,19)==1,:)=[];

nframes = size(data,1);
[C, ia, ic] = unique(data(:,1),'first');
tracklength = diff([ia;size(data,1)]);

data_opt = zeros(nframes,23);

% create mask & gutter
xi = [10 10 1270 1270];
yi = [1014 10 10 1014];
mask=poly2mask(xi,yi,double(ypix),double(xpix));
mask = uint8(mask);
gutter = find(mask==0);

checknum = 50;
nonlintic=tic;
% for iframe=1:nframes
    iframe=1;
    frametime = 0;
    for trackid=1:size(tracklength,1);
        
        new = 1;
        check = 1;

        for itrack=1:tracklength(trackid);
            
            flag = zeros(1,4);

            frame = tic;
            display(sprintf('Current track: %d out of %d tracks',trackid,size(tracklength,1)));
            display(sprintf('Current frame: %d out of %d in total',iframe,nframes));
            display(sprintf('Current frame: %d out of %d in track',itrack,tracklength(trackid)));

            ncams = size(cpvs,2);
            clusters = struct('ind',{[],[],[],[]},'indmod',{[],[],[],[]});

            posinfile=data(iframe,[8,11,14,17]);
            center_2d = struct('center',{data(iframe,6:7),data(iframe,9:10),data(iframe,12:13),data(iframe,15:16)});
            for icam=1:ncams
                nclusters=1;
                cpvfile = cpvs(icam);
                frewind(cpvfile);
                [cpvheader,rcnt] = fread(cpvfile,20,'*uint8');
                if rcnt ~= 20
                    warning('Error reading cpv file header');
                end

                xpix = uint32(cpvheader(5))+uint32(cpvheader(6))*256;
                ypix = uint32(cpvheader(7))+uint32(cpvheader(8))*256;
                fseek(cpvfile,posinfile(icam),'bof');
                currentFrameNum = fread(cpvfile,1,'*uint32');
    %             if currentFrameNum ~= cframe(icam)
    %                 error('frame number does not match');
    %             end
                [p, rcnt] = fread(cpvfile,4,'uint8=>uint32');
                pixcnt = bitshift(p(4),14)+bitshift(p(3),6)+bitshift(p(2),-2);
                pixlist = zeros(pixcnt,3,'uint32');
                dat = fread(cpvfile,[4,pixcnt],'uint8=>uint32');
                last3bitarray=uint32(07)+zeros(1,pixcnt,'uint32');
                pixlist(:,3)= dat(1,:);
                %pixel locations range [0,xpix-1] by [0,ypix -1]
                pixlist(:,1) = dat(2,:)+bitshift(bitand(dat(3,:),last3bitarray),8); % decoding x position
                pixlist(:,2) = bitshift(dat(3,:),-3) + bitshift(dat(4,:),5);% decoding y position
                img = zeros(ypix,xpix,'uint8');
                for i=1:pixcnt
                    if(pixlist(i,1) < xpix && pixlist(i,1) > 0 && pixlist(i,2) < ypix && pixlist(i,2) > 0)
                        img(pixlist(i,2)+1,pixlist(i,1)+1)=pixlist(i,3);
                    end
                end
                imgb = im2bw(img,0.01); % .25~64 .20~52 .175~43 filter out particles
%                 img = gpuArray(img);
%                 imgb = gpuArray(imgb);
                pix_info=regionprops(imgb,img,'Centroid','Area','PixelValues','PixelList','PixelIdxList');
                center=cat(1,pix_info.Centroid);
                area=cat(1,pix_info.Area);
                ind = find(area>cls_min(icam) & area<cls_max(icam));
                if ~isempty(ind)
                    for idz=size(ind,1):-1:1
                        if (~isempty(intersect(pix_info(ind(idz)).PixelIdxList,gutter)));
                            ind(idz)=[];
                        end
                    end
                end
                if ~isempty(ind)
                    for idz=1:size(ind,1)
                        clusters(icam).ind(nclusters:nclusters+size(pix_info(ind(idz)).PixelList,1)-1,1:2)=pix_info(ind(idz)).PixelList; 
                        clusters(icam).ind(nclusters:nclusters+size(pix_info(ind(idz)).PixelList,1)-1,3)=pix_info(ind(idz)).PixelValues;
                        nclusters = nclusters+size(pix_info(ind(idz)).PixelList,1);
                    end
                end
                if size(clusters(icam).ind,1)<cls_treshold
                    flag(icam) = 1;
                end
                clusters(icam).indmod=clusters(icam).ind;
            end
            % 0 if required on all cameras
            if size(flag(flag>0),2)>0
                data_opt(iframe,:) = [data(iframe,1) size(clusters(1).ind,1) size(clusters(2).ind,1) size(clusters(3).ind,1) size(clusters(4).ind,1) cls_treshold 0 data(iframe,5:end) -1];
                iframe=iframe+1;
                continue;
            end
            
            % Now we do a nonlinear least squares fit to find the orientation
            nonlin0=tic;
                        
            if new
                eul = rand(1,3)*2*pi;
                cntr = data(iframe,2:4);
                eul_cntr = [eul cntr];
            end
            
            if ~mod(itrack,checknum)
                check = 1;
            end
            
            [eul_cntr,fval,out,xflag] = sk_nonlinearopt_gpu(eul_cntr,clusters,camParaCalib,model,ICS_treshold,check);

            frametime=frametime+toc(frame);
            display(sprintf('\tThis frame: %f \n\tAverage time/frame: %f sec \n\tExpected to finish: %f min',toc(nonlin0), frametime/iframe,((nframes-iframe)*(frametime/iframe))/60));
            data_opt(iframe,:) = [data(iframe,1) eul_cntr(4:6) eul_cntr(1:3) data(iframe,5:end) fval];

            iframe=iframe+1;
            new = 0;
            check = 0;
        end
    %     save([savepath filename '_opt_temp.mat'], 'data_opt');
    end
% end

nonlint=toc(nonlintic);
display(sprintf('Nonlinear fit total time: %f seconds',nonlint));

save([savepath filename '_opt.mat'], 'data_opt');
% write2_gdf(data_opt,[savepath filename '_opt.gdf']);
fclose('all');