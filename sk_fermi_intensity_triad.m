function int_out=sk_fermi_intensity_triad(pix_xy,points_proj,I0)
% Inputs:  pix_xy:  (N_pixels by 2) matrix containing 2D (x,y) coordinates of bright pixels from data images
%     points_proj:  (N_arms*3 by 2) matrix containing 2D (x,y) coordinates of beginning, end and radius at midpoint of each arm
%                    points_proj(1:3,:) arms (endpoints)
%                    points_proj(4:6,:) cntr (startingpoints, center of triads)
%                    points_proj(7:9,:) rad (radius at midpoint of each arm)
%              I0:   Maximum intensity that should be used in returned intensity
% Outputs:int_out:  (N_pixels by 1) matrix with the intensity at each pixel in pix_xy

narms = 3; % number of arms

ends = points_proj(1:3,:); % endpoints of each arm
begs = points_proj(4:6,:); % beginning points of each arm (particle center)
rads = points_proj(7:9,:); % radius of each arm

npix = size(pix_xy,1);
intensity = zeros(npix,narms,'single');

mid = 0.5*(ends+begs); % center point of each arm
arml = sqrt(sum((ends-begs).^2,2));  % length of each arm in pixel
armw = sqrt(sum((rads-mid).^2,2)); % radius of each arm in pixel
theta = atan2(begs(:,2)-ends(:,2),begs(:,1)-ends(:,1));  % angle of each arm

for i=1:narms
    rot2D = [cos(theta(i)),sin(theta(i));-sin(theta(i)),cos(theta(i))];
    arm_co = abs((rot2D * (pix_xy-repmat(mid(i,:),npix,1))')');
    intensity(:,i) = I0*(1./(1+exp((arm_co(:,2)-armw(i))./(0.5*armw(i))))).*(1./(1+exp((arm_co(:,1)-.5*arml(i))./(0.025*arml(i)))));

%     [TX,TY] = meshgrid(-2*arml(i):2*arml(i), -2*armw(i):2*armw(i));
%     TI = double (I0*(1./(1+exp((abs(TY))-armw(i)./(0.5*armw(i))))).*(1./(1+exp((abs(TX)-.5*arml(i))./(0.025*arml(i))))));
%     figure(69);
%     surf(TX,TY,TI);
%     axis equal;
%     view(2);

end

int_out = max(intensity,[],2);

end
