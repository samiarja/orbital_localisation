%% Visualise the loss function output and plot the warped event for each selected region
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

VelArray = [-10:.01:10];
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
stride          = timeWindow/2;
stride_array    = 1:stride:e.t(end);

for iVelx = 1:nVel
    vx = VelArray(iVelx);
    for x = 1:xs
        xSum(x) = sum(round(e.x+vx*e.t/1e5) == x);
    end
    xStd(iVelx) = std(xSum);
end

for iVely = 1:nVel
    vy = VelArray(iVely);
    for y = 1:ys
        ySum(y) = sum(round(e.y+vy*e.t/1e5) == y);
    end
    yStd(iVely) = std(ySum);
end

[Vx,indx] = max(xStd);
[Vy,indy] = max(yStd);

velocity_x = VelArray(indx);
velocity_y = VelArray(indy);

VelX = [VelX;velocity_x];
VelY = [VelY;velocity_y];

figure(678789);
subtightplot(4,5,[1 16])
plot_shaded(VelArray(~isnan(yStd)),yStd(~isnan(yStd)));title("Vy");
xlabel("Velocity range");ylabel("Std");hold on;
V1 = -8;
plot(V1,yStd(find(VelArray==V1)),"og","MarkerSize",5);
plot(V1,yStd(find(VelArray==V1)),"og","MarkerSize",10);
plot(V1,yStd(find(VelArray==V1)),"og","MarkerSize",15);grid on
text(V1,yStd(find(VelArray==V1)),"\leftarrow" + num2str(V1),'FontSize',14)
plot(VelArray(indy),yStd(indy),"or","MarkerSize",5);
plot(VelArray(indy),yStd(indy),"or","MarkerSize",10);
plot(VelArray(indy),yStd(indy),"or","MarkerSize",15);grid on
text(VelArray(indy),yStd(indy),"\leftarrow" + num2str(VelArray(indy)),'FontSize',14)
V2 = -2;
plot(V2,yStd(find(VelArray==V2)),"ob","MarkerSize",5);
plot(V2,yStd(find(VelArray==V2)),"ob","MarkerSize",10);
plot(V2,yStd(find(VelArray==V2)),"ob","MarkerSize",15);grid on
text(V2,yStd(find(VelArray==V2)),"\leftarrow" + num2str(V2),'FontSize',14)
V3 = 5;
plot(V3,yStd(find(VelArray==V3)),"om","MarkerSize",5);
plot(V3,yStd(find(VelArray==V3)),"om","MarkerSize",10);
plot(V3,yStd(find(VelArray==V3)),"om","MarkerSize",15);grid on
text(V3,yStd(find(VelArray==V3)),"\leftarrow" + num2str(V3),'FontSize',14)
subtightplot(4,5,2)
scatter(e.x,round(e.y+V1*e.t/1e5),'.g','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,7)
scatter(e.x,round(e.y+V2*e.t/1e5),'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,12)
scatter(e.x,round(e.y+VelArray(indy)*e.t/1e5),'.r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,17)
scatter(e.x,round(e.y+V3*e.t/1e5),'.m','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

subtightplot(4,5,[3 18])
plot_shaded(VelArray(~isnan(xStd)),xStd(~isnan(xStd)));title("Vx");
xlabel("Velocity range");ylabel("Std");hold on;
V1 = -8;
plot(V1,xStd(find(VelArray==V1)),"og","MarkerSize",5);
plot(V1,xStd(find(VelArray==V1)),"og","MarkerSize",10);
plot(V1,xStd(find(VelArray==V1)),"og","MarkerSize",15);grid on
text(V1,xStd(find(VelArray==V1)),"\leftarrow" + num2str(V1),'FontSize',14)
plot(VelArray(indx),xStd(indx),"or","MarkerSize",5);
plot(VelArray(indx),xStd(indx),"or","MarkerSize",10);
plot(VelArray(indx),xStd(indx),"or","MarkerSize",15);grid on
text(VelArray(indx),xStd(indx),"\leftarrow" + num2str(VelArray(indx)),'FontSize',14)
V2 = -3;
plot(V2,xStd(find(VelArray==V2)),"ob","MarkerSize",5);
plot(V2,xStd(find(VelArray==V2)),"ob","MarkerSize",10);
plot(V2,xStd(find(VelArray==V2)),"ob","MarkerSize",15);grid on
text(V2,xStd(find(VelArray==V2)),"\leftarrow" + num2str(V2),'FontSize',14)
V3 = 5;
plot(V3,xStd(find(VelArray==V3)),"om","MarkerSize",5);
plot(V3,xStd(find(VelArray==V3)),"om","MarkerSize",10);
plot(V3,xStd(find(VelArray==V3)),"om","MarkerSize",15);grid on
text(V3,xStd(find(VelArray==V3)),"\leftarrow" + num2str(V3),'FontSize',14)
subtightplot(4,5,4)
scatter(e.y,round(e.x+V1*e.t/1e5),'.g','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,9)
scatter(e.y,round(e.x+V2*e.t/1e5),'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,14)
scatter(e.y,round(e.x+VelArray(indx)*e.t/1e5),'.r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,19)
scatter(e.y,round(e.x+V3*e.t/1e5),'.m','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

e.warpedx  = round(e.x+VelArray(indx)*e.t/1e5);
e.warpedy  = round(e.y+VelArray(indy)*e.t/1e5);
subtightplot(4,5,[5 10])
scatter(e.x,e.y,'.r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
subtightplot(4,5,[15 20])
scatter(e.warpedx,e.warpedy,'.r','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

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

figure(686);
subtightplot(1,2,1)
scatter3(e.x,e.y,e.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
subtightplot(1,2,2)
scatter(e.x,e.y,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

disp("Dataset Generated");

%%%%%%%%%%%%%%%%%%%% Event Warping algorithm %%%%%%%%%%%%%%%%%%
disp("Start Event Warping algorithm");
stride_counter = 0;
finalSTD = nan(numel(1:10:1000),1);
xFinalStd = nan(1,numel(finalSTD));
xFinalSum = nan(1,xs);
yFinalStd = nan(1,numel(finalSTD));
yFinalSum = nan(1,ys);

for stdidx = 10
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
        counter
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

%% Apply rotation on the events
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


VelArray = [-10:.01:10];
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
stride          = timeWindow/2;
stride_array    = 1:stride:e.t(end);

counterx = 0;
angleRange = 0.001:0.001:1;
finalSTDx = nan(nVel, numel(0.001:0.001:1));

for angle = 0.001:0.001:1
    counterx = counterx+1;
    counterx
    theta = angle;
    Rotation_mat = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    rot_coor = [e.x; e.y];
    rotxy = Rotation_mat*rot_coor;
    e.x = rotxy(1,:);e.y = rotxy(2,:);
    
    for iVelx = 1:nVel
        vx = VelArray(iVelx);
        for x = 1:xs
            xSum(x) = sum(round(e.x+vx*e.t/1e5) == x);
        end
%         xStd(iVelx) = std(xSum);
        finalSTDx(iVelx,counterx) = std(xSum);
%         figure(244347767);
%         scatter(e.y,round(e.x+vx*e.t/1e5),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
%         xlabel("X");ylabel("Y");
    end
end

counter = 0;
angleRange = 0.001:0.001:1;
finalSTD = nan(nVel, numel(0.001:0.001:1));
for angle = 0.001:0.001:1
    counter = counter+1;
    counter
    theta = angle;
    Rotation_mat = [cos(theta) -sin(theta);sin(theta) cos(theta)];
    rot_coor = [e.x; e.y];
    rotxy = Rotation_mat*rot_coor;
    e.x = rotxy(1,:);e.y = rotxy(2,:);
    
    for iVely = 1:nVel
        vy = VelArray(iVely);
        for y = 1:ys
            ySum(y) = sum(round(e.y+vy*e.t/1e5) == y);
        end
        finalSTD(iVely,counter) = std(ySum);
%         figure(56453534);
%         scatter(e.x,round(e.y+vy*e.t/1e5),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
    end
end

vx=-0.0300;e.x = round(e.x+vx*e.t/1e5);
vy=0.0300;e.y = round(e.y+vy*e.t/1e5);

theta = 0.05;
Rotation_mat = [cos(theta) -sin(theta);sin(theta) cos(theta)];
rot_coor = [e.x; e.y];
rotxy = Rotation_mat*rot_coor;
e.x = rotxy(1,:);e.y = rotxy(2,:);

figure(6779);
scatter(e.x,e.y,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("Vy")

figure(33545);
scatter(e.y,e.x,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("Vx")

figure(46565);
scatter3(e.x,e.y,e.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);

figure(879789789);
subtightplot(2,3,1)
scatter(e.y,e.x,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("Unwarped Events")
subtightplot(2,3,2)
scatter(ewarped.y,ewarped.x,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("Warped Events")
subtightplot(2,3,4)
maximum = nan(numel(finalSTD(1,:)),1);
for idx = 1:numel(finalSTD(1,:))
   maximum(idx,1) = max(finalSTD(:,idx));
end
[value,index] = max(maximum);
plot(finalSTD);grid on;xlim([0 400]);hold on
plot(finalSTD(:,index),"LineWidth",2);
plot(150,value,"or","MarkerSize",5);
plot(150,value,"or","MarkerSize",10);
plot(150,value,"or","MarkerSize",15);grid on
text(150,value,"\leftarrow" + "Angle_{x}: " + num2str(angleRange(index)),'FontSize',14)
subtightplot(2,3,5)
maximum = nan(numel(finalSTDx(1,:)),1);
for idx = 1:numel(finalSTDx(1,:))
   maximum(idx,1) = max(finalSTDx(:,idx));
end
[value,index] = max(maximum);
plot(finalSTDx);grid on;xlim([0 400]);hold on
plot(finalSTDx(:,index),"LineWidth",2);
plot(147,value,"or","MarkerSize",5);
plot(147,value,"or","MarkerSize",10);
plot(147,value,"or","MarkerSize",15);grid on
text(147,value,"\leftarrow" + "Angle_{y}: " + num2str(angleRange(index)),'FontSize',14)
subtightplot(2,3,[3 6])
scatter3(ewarped.x,ewarped.y,ewarped.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("3D Warped Events");

%% visualise window size with respect to window and stride size
load('test_files/window_wrt_std.mat')

figure(6786);
subplot(2,6,1)
scatter(e_final.warped{1, 1}.warpedy,e_final.warped{1, 1}.warpedx,'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("X (px)");
subplot(2,6,2)
scatter(e_final.warped{1, 2}.warpedy,e_final.warped{1, 2}.warpedx,'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("X (px)");
subplot(2,6,3)
scatter(e_final.warped{1, 3}.warpedy,e_final.warped{1, 3}.warpedx,'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("X (px)");
subplot(2,6,4)
scatter(e_final.warped{1, 4}.warpedy,e_final.warped{1, 4}.warpedx,'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("X (px)");
subplot(2,6,5)
scatter(e_final.warped{1, 10}.warpedy,e_final.warped{1, 10}.warpedx,'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("X (px)");
subplot(2,6,6)
scatter(e_final.warped{1, end}.warpedy,e_final.warped{1, end}.warpedx,'.b','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);
xlabel("X (px)");ylabel("X (px)");
subplot(2,6,[7 12])
plot_shaded(1:numel(yFinalStd(~isnan(yFinalStd))),yFinalStd(~isnan(yFinalStd)));xlabel("Stride");ylabel("Std");
% subplot(2,6,[10 12])
% plot_shaded(1:numel(yFinalStd(~isnan(yFinalStd))),yFinalStd(~isnan(yFinalStd)));

%% investigate the speed in relation to std
%%%%%%%%%%%%%%%%%%%%%% GENERATE THE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SPEED = 1;
ys = 500;
xs = 500;
diameter = 3;
nTime = 4e6;
nBackEvent = 4000;
nStarEvent = 4000;
lambda0 = nBackEvent/nTime;
lambda1 = nStarEvent/nTime;

backSpike = rand(1,nTime) < lambda0;
starSpike = rand(1,nTime) < lambda1;
centrex_start = 50;centrex_end = 50;
centrey_start = 50;centrey_end = 75;
vx_gt = (centrex_end - centrex_start) / nTime;
vy_gt = (centrey_end - centrey_start) / nTime;
speedC = 0;

xStdF = nan(1,numel(0:0.00001:0.00009));
yStdF = nan(1,numel(0:0.00001:0.00009));
e_final = [];
for S = 0:0.00001:0.00009
    speedC = speedC + 1;
    S
    kx = S;%[0 9e-5]; %9e-5;
    ky = S;%[0 9e-5]; %9e-5;
    
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
%     figure(686);
%     scatter3(e.x,e.y,e.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
    e_final.input{speedC}  = [e.x; e.y ;e.t];
    
    %%%%%%%%%%%%%%%%%%%%%% WARP THE EVENTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    VelArray = [-2:.01:2];
    nVel = numel(VelArray);
    xStd = nan(1,nVel);
    xSum = nan(1,xs);
    yStd = nan(1,nVel);
    ySum = nan(1,ys);
    
    
    xSumF = nan(1,xs);
    ySumF = nan(1,ys);
    
    counter             = 0;
    VelX                = [];
    VelY                = [];
    warped_events_x     = [];
    warped_events_y     = [];
    warpedT_events_x    = [];
    warpedT_events_y    = [];
    oldMat              = [];
    timeWindow          = 100000; % 1e6;
    stride              = timeWindow/5;
    stride_array        = 1:stride:e.t(end);
    xSumWarped          = nan(1,xs);
    ySumWarped          = nan(1,ys);
    coor_center         = nan(numel(stride_array),3);
    
    for ttd = stride_array-1
        counter = counter + 1;
        ii = find(e.t> stride_array(counter) & e.t<(timeWindow+stride_array(counter)));
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
                ii(overlapIdx) = [];
                
                warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e5)'];
                warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e5)'];

                for x = 1:xs
                    xSumWarped(x) = sum(round(e.x(ii)+velocity_x*e.t(ii)/1e5) == x);
                end
                for y = 1:ys
                    ySumWarped(y) = sum(round(e.y(ii)+velocity_y*e.t(ii)/1e5) == y);
                end
                [valx, indX] = max(xSumWarped);
                [valt, indY] = max(ySumWarped);
                coor_center(counter-1,1) = indX;coor_center(counter-1,2) = indY;coor_center(counter-1,3) = counter-1;

                diffX = coor_center(1,1) - coor_center(counter-1,1);
                diffY = coor_center(1,2) - coor_center(counter-1,2);

                warpedT_events_x = [warpedT_events_x;round(e.x(ii)+velocity_x*e.t(ii)/1e5)'+diffX];
                warpedT_events_y = [warpedT_events_y;round(e.y(ii)+velocity_y*e.t(ii)/1e5)'+diffY];
                
            end
        end
    end
    
    e.warpedx  = warped_events_x;
    e.warpedy  = warped_events_y;
    
    e.warpedTx  = warpedT_events_x;
    e.warpedTy  = warpedT_events_y;
    
    for x = 1:xs
        xSumF(x) = sum(e.warpedTx == x);
    end
    xStdF(speedC) = std(xSumF);
    for y = 1:ys
        ySumF(y) = sum(e.warpedTy == y);
    end
    yStdF(speedC) = std(ySumF);
    
    e_final.output{speedC} = [e.warpedTx e.warpedTy];
    
end

figure(346456);
subplot(3,5,1)
scatter3(e_final.input{1, 1}(2,:),e_final.input{1, 1}(1,:),e_final.input{1, 1}(3,:),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("No speed")
subplot(3,5,2)
scatter3(e_final.input{1, 3}(2,:),e_final.input{1, 3}(1,:),e_final.input{1, 3}(3,:),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("Speed: 1e-4")
subplot(3,5,3)
scatter3(e_final.input{1, 5}(2,:),e_final.input{1, 5}(1,:),e_final.input{1, 5}(3,:),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("Speed: 5e-4")
subplot(3,5,4)
scatter3(e_final.input{1, 7}(2,:),e_final.input{1, 7}(1,:),e_final.input{1, 7}(3,:),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("Speed: 7e-4")
subplot(3,5,5)
scatter3(e_final.input{1, 10}(2,:),e_final.input{1, 10}(1,:),e_final.input{1, 10}(3,:),'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2)
xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("Speed: 9e-4")
subplot(3,5,[6 10])
plot_shaded(1:numel(xStdF(~isnan(xStdF))),xStdF(~isnan(xStdF))); hold on
plot_shaded(1:numel(yStdF),yStdF);xlabel("Speed");ylabel("Std");grid on;hold on;
plot(1,xStdF(1),"or","MarkerSize",5);
plot(1,xStdF(1),"or","MarkerSize",10);
plot(1,xStdF(1),"or","MarkerSize",15);grid on
plot(3,xStdF(3),"or","MarkerSize",5);
plot(3,xStdF(3),"or","MarkerSize",10);
plot(3,xStdF(3),"or","MarkerSize",15);grid on
plot(5,xStdF(5),"or","MarkerSize",5);
plot(5,xStdF(5),"or","MarkerSize",10);
plot(5,xStdF(5),"or","MarkerSize",15);grid on
plot(7,xStdF(7),"or","MarkerSize",5);
plot(7,xStdF(7),"or","MarkerSize",10);
plot(7,xStdF(7),"or","MarkerSize",15);grid on
plot(10,xStdF(10),"or","MarkerSize",5);
plot(10,xStdF(10),"or","MarkerSize",10);
plot(10,xStdF(10),"or","MarkerSize",15);grid on
subtightplot(3,5,11)
imagesc(accumulate_render(e_final.output{1, 1}(:,1),e_final.output{1, 1}(:,2)).^(1 / 7));colormap(magma(100));colorbar;axis off
subtightplot(3,5,12)
imagesc(accumulate_render(e_final.output{1, 3}(:,1),e_final.output{1, 3}(:,2)).^(1 / 7));colormap(magma(100));colorbar;axis off
subtightplot(3,5,13)
imagesc(accumulate_render(e_final.output{1, 5}(:,1),e_final.output{1, 5}(:,2)).^(1 / 7));colormap(magma(100));colorbar;axis off
subtightplot(3,5,14)
imagesc(accumulate_render(e_final.output{1, 7}(:,1),e_final.output{1, 7}(:,2)).^(1 / 7));colormap(magma(100));colorbar;axis off
subtightplot(3,5,15)
imagesc(accumulate_render(e_final.output{1, 10}(:,1),e_final.output{1, 10}(:,2)).^(1 / 7));colormap(magma(100));colorbar;axis off

% figure(76767);
% subplot(1,4,1)
% scatter3(e.warpedx, e.warpedy,e.t,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
% xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("3D Warp - Non Motion Compensated")
% subplot(1,4,2)
% scatter(e.warpedx, e.warpedy,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
% xlabel("X (px)");ylabel("Y (px)");title("2D Warp - Non Motion Compensated")
% subplot(1,4,3)
% scatter3(e.warpedTx, e.warpedTy,e.t(1:numel(e.warpedTy)),'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
% xlabel("X (px)");ylabel("Y (px)");zlabel("Time \mu s");title("3D Warp - Motion Compensated")
% subplot(1,4,4)
% scatter(e.warpedTx, e.warpedTy,'.','MarkerFaceAlpha',.5,'MarkerEdgeAlpha',.5);hold on
% xlabel("X (px)");ylabel("Y (px)");title("2D Warp - Motion Compensated")


