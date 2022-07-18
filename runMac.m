%% apply window and stride and investigate the variation with different window/stride size
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

e = []; idx = 0;e_final = [];
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
% plot3(e.x,e.y,e.t,'.')

disp("Dataset Generated");

%%%%%%%%%%%%%%%%%%%% Event Warping algorithm %%%%%%%%%%%%%%%%%%
disp("Start Event Warping algorithm");
stride_counter = 0;
finalSTD = nan(numel(1:10:1000),1);
xFinalStd = nan(1,numel(finalSTD));
xFinalSum = nan(1,xs);
yFinalStd = nan(1,numel(finalSTD));
yFinalSum = nan(1,ys);

for stdidx = 110:10:1000
    stdidx
    stride_counter = stride_counter + 1;
    
    VelArray = [-2:.01:2];
    nVel = numel(VelArray);
    xStd = nan(1,nVel);
    xSum = nan(1,xs);
    yStd = nan(1,nVel);
    ySum = nan(1,ys);
    
    counter         = 10;
    VelX            = [];
    VelY            = [];
    warped_events_x = [];
    warped_events_y = [];
    oldMat          = [];
    timeWindow      = 100000; % 1e6;
    
    stride          = timeWindow/stdidx;
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
    
    e_final.warped{stride_counter} = e;
    
    for x = 1:xs
        xFinalSum(x) = sum(e.warpedx == x);
    end
    xFinalStd(stride_counter) = std(xFinalSum);
    
    for y = 1:ys
        yFinalSum(y) = sum(e.warpedy == y);
    end
    yFinalStd(stride_counter) = std(yFinalSum);
end

figure(68678);
subplot(1,2,1)
plot(xFinalStd,'-b');grid on;hold on;
subplot(1,2,2)
plot(yFinalStd,'-r');grid on