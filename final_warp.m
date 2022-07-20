%%%%%%%%%%%%%%%%%%%%%% GENERATE THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPEED = 1;
ys = 100;
xs = 100;
diameter = 3;
nTime = 4e6;
nBackEvent = 2000;
nStarEvent = 4000;
lambda0 = nBackEvent/nTime;
lambda1 = nStarEvent/nTime;

backSpike = rand(1,nTime) < lambda0;
starSpike = rand(1,nTime) < lambda1;
centrex_start = 50;centrex_end = 50;
centrey_start = 50;centrey_end = 75;
vx_gt = (centrex_end - centrex_start) / nTime;
vy_gt = (centrey_end - centrey_start) / nTime;
kx = 2e-5;
ky = 5e-6;

Vx_init = vx_gt;Vy_init = vy_gt;
e = []; idx = 0;
for t = 1:nTime
    if SPEED
        vx_gt = (centrex_end - centrex_start + kx*t) / nTime;
        vy_gt = (centrey_end - centrey_start + ky*t) / nTime;
    end
    starX = t*vx_gt + centrex_start;
    starY = t*vy_gt + centrey_start;
    if backSpike(t)
        idx = idx+1;
        x = randi(xs,1);
        y = randi(ys,1);
        e.x(idx) = round(x);
        e.y(idx) = round(y);
        e.p(idx) = 1;
        e.t(idx) = t;
    end
    
    if starSpike(t)
        idx = idx+1;
        x = randi(diameter,1) + starX;
        y = randi(diameter,1) + starY;
        e.x(idx) = round(x);
        e.y(idx) = round(y);
        e.p(idx) = 1;
        e.t(idx) = t;
    end
end
Vx_final = vx_gt;Vy_final = vy_gt;
% figure(686);
% scatter3(e.x,e.y,e.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)

%%%%%%%%%%%%%%%%%%%%%% WARP THE EVENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VelArray = [-2:.01:2];
nVel = numel(VelArray);
xStd = nan(1,nVel);
xSum = nan(1,xs);
yStd = nan(1,nVel);
ySum = nan(1,ys);

counter             = 0;
VelX                = [];
VelY                = [];
warped_events_x     = [];
warped_events_y     = [];
warpedT_events_x    = [];
warpedT_events_y    = [];
oldMat              = [];
e_final             = [];
timeWindow          = 100000; % 1e6;
stride              = timeWindow/5;
stride_array        = 1:stride:e.t(end);
xSumWarped          = nan(1,xs);
ySumWarped          = nan(1,ys);
coor_center         = nan(numel(stride_array),3);

for ttd = stride_array-1
    counter = counter + 1;
    counter
    ii = find(e.t> stride_array(counter) & e.t<(timeWindow+stride_array(counter)));
    if ~isempty(ii)
        for iVelx = 1:nVel
            vx = VelArray(iVelx);
            for x = 1:xs
                xSum(x) = summation(round(e.x(ii)+vx*e.t(ii)/1e5) == x);
            end
            xStd(iVelx) = std(xSum);
        end
        
        for iVely = 1:nVel
            vy = VelArray(iVely);
            for y = 1:ys
                ySum(y) = summation(round(e.y(ii)+vy*e.t(ii)/1e5) == y);
            end
            yStd(iVely) = std(ySum);
        end
        
        [Vx,indx] = max(xStd);
        [Vy,indy] = max(yStd);
        
        velocity_x = VelArray(indx);
        velocity_y = VelArray(indy);
        VelX = [VelX;velocity_x];
        VelY = [VelY;velocity_y];
        
        if counter == 1
            warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e5)'];
            warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e5)'];
        else
            oldMat = find(e.t> stride_array(counter-1) & e.t<(timeWindow+stride_array(counter-1)));
            [overlapVal,overlapIdx] = intersect(ii,oldMat);
            ii(overlapIdx) = [];
            
            e.warpedx  = round(e.x(ii)+velocity_x*e.t(ii)/1e5)';
            e.warpedy  = round(e.y(ii)+velocity_y*e.t(ii)/1e5)';
            e.warpedt  = (e.t(ii)/1e5)';
            
            for x = 1:xs
                xSumWarped(x) = summation(e.warpedx == x);
            end
            for y = 1:ys
                ySumWarped(y) = summation(e.warpedy == y);
            end
            [valx, indX] = max(xSumWarped);
            [valt, indY] = max(ySumWarped);
            coor_center(counter-1,1) = indX;coor_center(counter-1,2) = indY;coor_center(counter-1,3) = counter-1;
            
            diffX = coor_center(1,1) - coor_center(counter-1,1);
            diffY = coor_center(1,2) - coor_center(counter-1,2);
            
            e.warpedTx = e.warpedx + diffX;
            e.warpedTy = e.warpedy + diffY;
            
            warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e5)'];
            warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e5)'];
            
            warpedT_events_x = [warpedT_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e5)'+diffX];
            warpedT_events_y = [warpedT_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e5)'+diffY];
            
%             figure(76767);
%             subplot(1,4,1)
%             scatter3(e.warpedx, e.warpedy,e.warpedt,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
%             xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("3D Warp - Non Motion Compensated")
%             subplot(1,4,2)
%             scatter(e.warpedx, e.warpedy,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
%             xlabel("X (px)");ylabel("Y (px)");title("2D Warp - Non Motion Compensated")
%             subplot(1,4,3)
%             scatter3(e.warpedTx, e.warpedTy,e.warpedt,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
%             xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("3D Warp - Motion Compensated")
%             subplot(1,4,4)
%             scatter(e.warpedTx, e.warpedTy,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
%             xlabel("X (px)");ylabel("Y (px)");title("2D Warp - Motion Compensated")
        end
    end
end
e.finalX = warped_events_x;
e.finalY = warped_events_y;

e.finalTX = warpedT_events_x;
e.finalTY = warpedT_events_y;

figure(76767);
subplot(1,4,1)
scatter3(e.finalX, e.finalY,1:numel(e.finalY),'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("3D Warp - Non Motion Compensated")
subplot(1,4,2)
scatter(e.finalX, e.finalY,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
xlabel("X (px)");ylabel("Y (px)");title("2D Warp - Non Motion Compensated")
subplot(1,4,3)
scatter3(e.finalTX, e.finalTY,1:numel(e.finalTY),'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("3D Warp - Motion Compensated")
subplot(1,4,4)
scatter(e.finalTX, e.finalTY,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
xlabel("X (px)");ylabel("Y (px)");title("2D Warp - Motion Compensated")

%%%%%%%%%%%%%%%%%%%%%%%%%%% ACCUMULATE AND RENDER %%%%%%%%%%%
sigma = 0.5;
padding = double(ceil(6.0 * sigma));
kernel_indices = zeros(padding * 2 + 1, padding * 2 + 1, 2);

for y = 1:padding * 2 + 1
    for x = 1:padding * 2 + 1
        kernel_indices(y, x, 1) = x - padding;
        kernel_indices(y, x, 2) = y - padding;
    end
end
x_minimum = min(e.finalTX);
y_minimum = min(e.finalTY);

xs = e.finalTX - x_minimum + padding;
ys = e.finalTY - y_minimum + padding;

pixels = zeros(ceil(max(ys)) + padding + 1,ceil(max(xs)) + padding + 1);

xis = round(xs);
yis = round(ys);

xfs = xs - xis;
yfs = ys - yis;

sigma_factor = -1.0 / (2.0 * sigma^2.0);

for i = 1:numel(xis)
    sumF = zeros(size(kernel_indices,1),size(kernel_indices,1),size(kernel_indices,3));
    for l = 1:size(kernel_indices,1)
        for j = 1:size(kernel_indices,1)
            summation = 0;
            for axis  = 1:size(kernel_indices,3)
                summation = summation + kernel_indices(l,j,axis);
                sumF(l,j,axis) = summation;
            end
        end
    end
    sumF = sumF.^ 2.0;
    finalSummation = (sumF(:,:,1) + sumF(:,:,2))*sigma_factor;
    pixels(yis(i)-padding+1:yis(i)+padding+1,xis(i)-padding+1:xis(i)+padding+1) = pixels(yis(i)-padding+1:yis(i)+padding+1,xis(i)-padding+1:xis(i)+padding+1)+...
        exp(finalSummation);
end

figure(78788);
imagesc(pixels.^(1 / 3));colormap(magma(100));colorbar



