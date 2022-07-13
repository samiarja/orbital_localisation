%% Event warping of ISS data
load("./data/mat/20220125_Brittany_211945_2022-01-25_21-21-18_NADIR.mat")
e = struct("x",double(events(:,2)),"y",double(events(:,3)),"p",double(events(:,4)),"ts",double(events(:,1)));
e = FilterTD(e, 2);
e.t = e.ts;
e.t = e.t - e.t(1);

VelArray = [-2:.01:2];
nVel = numel(VelArray);
xStd = nan(1,nVel);
xSum = nan(1,xs);
yStd = nan(1,nVel);
ySum = nan(1,ys);

counter         = 0;
VelX            = [];
VelY            = [];
warped_events_x = [];
warped_events_y = [];
oldMat          = [];
timeWindow      = 10000000; % 1e6;
stride          = timeWindow/2;
stride_array    = 1:stride:e.t(end);

for ttd = stride_array-1
    counter = counter + 1;
    counter
    ii = find(e.t> stride_array(counter) & e.t<(timeWindow+stride_array(counter))); % with stride
%     ii = find(e.t>stride_array(counter) & e.t<(stride_array(counter)+stride)); % without stride
%     ii = 1:numel(e.t); % all events
    if ~isempty(ii)
        for iVelx = 1:nVel
            vx = VelArray(iVelx);
            for x = 1:xs
                xSum(x) = sum(round(e.x(ii)+vx*e.t(ii)/1e6) == x);
            end
            xStd(iVelx) = std(xSum);
        end
        
        for iVely = 1:nVel
            vy = VelArray(iVely);
            for y = 1:ys
                ySum(y) = sum(round(e.y(ii)+vy*e.t(ii)/1e6) == y);
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
            warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e6)'];
            warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e6)'];
        else
            oldMat = find(e.t> stride_array(counter-1) & e.t<(timeWindow+stride_array(counter-1)));
            [overlapVal,overlapIdx] = intersect(ii,oldMat);
%             newMatrixTest = ii;
            ii(overlapIdx) = [];
            warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e6)'];
            warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e6)'];
        end
    end
end
e.warpedx  = warped_events_x;
e.warpedy  = warped_events_y;

figure(45667);
subplot(3,3,[1 4])
scatter3(e.x,e.y,e.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("Input Events");
xlabel("X (px)");ylabel("Y (px)");zlabel("Time (\mu s)")
subplot(3,3,[2 5])
scatter3(e.warpedx,e.warpedy,e.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("Y (px)");zlabel("Time (\mu s)")
title("3D Warped Events");
subplot(3,3,3)
scatter(e.x,e.y,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("2D Raw Events (No Warp)");
xlabel("X (px)");ylabel("Y (px)");
subplot(3,3,6)
scatter(e.warpedx,e.warpedy,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("2D Warped Events");
xlabel("X (px)");ylabel("Y (px)");
subplot(3,3,[7 8])
plot_shaded(1:numel(VelX),VelX);hold on;
plot_shaded(1:numel(VelY),VelY);ylabel("Velocity (max of the std)");xlabel("Slice of event");
y1 = yline(-vx_gt*1e5,'-.','Vx\_gt','LineWidth',2);y1.LabelHorizontalAlignment = 'left';
y2 = yline(-vy_gt*1e5,'-.','Vy\_gt','LineWidth',2);y2.LabelHorizontalAlignment = 'left';
grid on;title("Motion Estimated Velocity \color{blue}Vx, \color{red}Vy");
subplot(3,3,9)
pts = linspace(min(e.warpedx), max(e.warpedx), 80);
N = histcounts2(e.warpedx, e.warpedy, pts, pts);
imagesc(pts, pts, N);colorbar;title("Warped Events HeatMap");colormap("gray");xlabel("X (px)");ylabel("Y (px)");
