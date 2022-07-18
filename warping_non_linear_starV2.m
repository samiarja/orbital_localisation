%% Generate dataset (Fake star moving across the sky)
SPEED = 1;
ys = 100;
xs = 300;
diameter = 3;
nTime = 4e6; % 4 sec
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
kx = 2e-5; % small value (velocity is function of time)
ky = 5e-6; % small value

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

figure(686);
plot3(e.x,e.y,e.t,'.')
figure(7879);
plot(-[Vx_init Vx_final],'-');hold on;
plot(-[Vy_init Vy_final],'-');grid on;legend("Vx","Vy");

% nEvents = numel(e.x);
% figure(2); clf;
% S = zeros(xs,ys); T = S; P = T;%-inf;
% ss = imagesc(S); colorbar; axis image;

% tau = 25e3;
% displayFreq = 3e3; % in units of time
% nextTimeSample = displayFreq;
% for idx = 1:round(nEvents)
%     x = e.x(idx);
%     y = e.y(idx);
%     t = e.t(idx);
%     p = e.p(idx);
%     T(x,y) = t;
%     P(x,y) = p;
%     if t > nextTimeSample
%         nextTimeSample = nextTimeSample + displayFreq;
%         S = P.*exp((T-t)/tau);
%         set(ss,'CData',S)
%         drawnow
%     end
% end
%% Event warping
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
timeWindow      = 100000; % 1e6;
stride          = timeWindow/50;
stride_array    = 1:stride:e.t(end);

for ttd = stride_array-1
    counter = counter + 1;
    ii = find(e.t> stride_array(counter) & e.t<(timeWindow+stride_array(counter))); % with stride
%     ii = find(e.t>stride_array(counter) & e.t<(stride_array(counter)+stride)); % without stride
%     ii = 1:numel(e.t); % all events
    if ~isempty(ii)
        for iVelx = 1:nVel
            vx = VelArray(iVelx);
            for x = 1:xs
                xSum(x) = sum(round(e.x(ii)+vx*e.t(ii)/1e5) == x);
            end
            xStd(iVelx) = std(xSum);
%             figure(15345345)
%             plot(e.y(ii),round(e.x(ii)+vy*e.t(ii)/1e6) ,'.')
%             figure(35345345)
%             plot(xSum,1:xs)
%             grid on;
%             axis image
%             drawnow
        end
        
        for iVely = 1:nVel
            vy = VelArray(iVely);
            for y = 1:ys
                ySum(y) = sum(round(e.y(ii)+vy*e.t(ii)/1e5) == y);
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
%             newMatrixTest = ii;
            ii(overlapIdx) = [];
            warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e5)'];
            warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e5)'];
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
scatter3(e.warpedx,e.warpedy,e.t(1:numel(e.warpedy)),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
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
