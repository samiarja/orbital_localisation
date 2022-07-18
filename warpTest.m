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
scatter3(ewarped.x,ewarped.y,ewarped.t,'.','MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2);title("3D Warped Events")
