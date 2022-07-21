%% Event warping with non linear velocity
%TODO: calculate Vx and Vy with std as a objective function
% calculate the new events based on the new velocity
% render the event

% load('data/mat/20220125_Brittany_211945_2022-01-25_21-21-18_NADIR.mat')
load('data/mat/FN034HRTEgyptB_NADIR.mat')
TD = struct("x",double(events(:,2)),"y",double(events(:,3)),"p",double(events(:,4)),"ts",double(events(:,1)));
TD = FilterTD(TD, 10);
TD.ts = TD.ts - TD.ts(1);
velocity_array_final_x = [];
velocity_array_final_y = [];
counter = 0;
dimensionX = 0;
dimensionY = 0;
warp = [];
warped_events = [];
pano_height = 1024;
pano_width = 2 * pano_height;
S = zeros(pano_height,pano_width);
time = 1e6:1e6:TD.ts(end);
timestamp = [];
new_x = [];
new_y = [];

for slice = time
    counter = counter + 1;
    tt = TD.ts;
    xx = TD.x;
    yy = TD.y;
    pp = TD.p;
    %     t0 = tt(1)+slice;
    t0 = time(counter)+slice;
    iFirst = find(tt<t0,1,'last');
    iLast = find(tt>(t0+1e6),1,'first');
    ii = iFirst:iLast;
    disp((double(tt(iLast) - tt(iFirst)))/1e6)
    tt(ii);
    
    if ~isempty(ii)
        nX = 240; nY = 180;
        
        % xSum = nan(1,nX);
        % for x = 1:nX
        %     xSum(x) = sum(round((xx(ii)+vx*tt(ii)/1e6)) == x);
        % end
        % figure(25345345)
        % plot(xSum)
        % grid on;
        %xStd(iVelx) = std(xSum)
        
        VelArray = [-2:.01:2];
        nVel = numel(VelArray);
        xStd = nan(1,nVel);
        xSum = nan(1,nX);
        yStd = nan(1,nVel);
        ySum = nan(1,nY);
        for iVelx = 1:nVel
            %%% loop through Vx
            vx = VelArray(iVelx);
            for x = 1:nX
                xSum(x) = sum(round(xx(ii)+vx*tt(ii)/1e6) == x);
            end
            xStd(iVelx) = std(xSum);
            %%% loop through Vy
            vy = VelArray(iVelx);
            for y = 1:nY
                ySum(y) = sum(round(yy(ii)+vy*tt(ii)/1e6) == y);
            end
            yStd(iVelx) = std(ySum);
%             figure(15345345)
%             plot(xx(ii),round(yy(ii)+vy*tt(ii)/1e6) ,'.')
%             figure(35345345)
%             plot(ySum,1:nY)
%             grid on;
%             axis image
%             drawnow
        end
        
        [Vx,indx] = max(xStd);
        [Vy,indy] = max(yStd);
        
        %     figure(53468)
        %     subplot(2,1,1)
        %     plot(VelArray,xStd);hold on
        %     plot(VelArray(indx),Vx,'or','MarkerSize',10);title("Vx");
        %     subplot(2,1,2)
        %     plot(VelArray,yStd);hold on
        %     plot(VelArray(indy),Vy,'or','MarkerSize',10);title("Vy");
        
        velocity_x = VelArray(indx);
        velocity_y = VelArray(indy);
        
        warped_events.x  = round(xx(ii)+velocity_x*tt(ii)/1e6);
        warped_events.y  = round(yy(ii)+velocity_y*tt(ii)/1e6);
        warped_events.p  = pp(ii);
        warped_events.ts = tt(ii);
        warped_events.vx = repelem(velocity_x,numel(ii))';
        warped_events.vy = repelem(velocity_y,numel(ii))';
        
        warp.warp{counter} =  warped_events;
        
        dimensionX = dimensionX + max(warped_events.x);
        dimensionY = dimensionY + max(warped_events.y);
        velocity_array_final_x=[velocity_array_final_x;velocity_x];
        velocity_array_final_y=[velocity_array_final_y;velocity_y];
        
        figure(67867);
        scatter(warped_events.x,warped_events.y,'.');
        
    end
end

for idx = 1:numel(warp.warp)
    if isempty(warp.warp)
        warp.warp(idx) = [];
    end
end

figure(789789);
plot(velocity_array_final_x);hold on
plot(velocity_array_final_y);

% combine all warp arrays into one single event stream
motion_compensated_imagex  = [];
motion_compensated_imagey  = [];
motion_compensated_imagep  = [];
motion_compensated_imagets = [];
motion_compensated_imagevx = [];
motion_compensated_imagevy = [];
for idx = 1:numel(warp.warp)
    motion_compensated_imagex  = [motion_compensated_imagex;warp.warp{1, idx}.x];
    motion_compensated_imagey  = [motion_compensated_imagey;warp.warp{1, idx}.y];
    motion_compensated_imagep  = [motion_compensated_imagep;warp.warp{1, idx}.p];
    motion_compensated_imagets = [motion_compensated_imagets;warp.warp{1, idx}.ts];
    motion_compensated_imagevx = [motion_compensated_imagevx;warp.warp{1, idx}.vx];
    motion_compensated_imagevy = [motion_compensated_imagevy;warp.warp{1, idx}.vx];
%     motion_compensated_imagevx = [motion_compensated_imagevx;repelem(warp.warp{1, idx}.vx,numel(warp.warp{1, idx}.x))'];
%     motion_compensated_imagevy = [motion_compensated_imagevy;repelem(warp.warp{1, idx}.vy,numel(warp.warp{1, idx}.x))'];
end
motion_compensated_image = struct("x",motion_compensated_imagex, ...
                                  "y",motion_compensated_imagey, ...
                                  "p",motion_compensated_imagep, ...
                                  "ts",motion_compensated_imagets, ...
                                  "vx",motion_compensated_imagevx, ...
                                  "vy",motion_compensated_imagevy);
%% visualisation
load('test_files/warp_nofilter_egypt.mat');
png_satellite_night=imageDatastore("data/Egypt/png_satellite/*.png");
velocityX = [];velocityY = [];
frame_counter = 0;
sigma = 0.7;
pts = linspace(-10, 100, 80);
panoramic_satellite = stitch("data/Egypt/png_event/");

% writerObj = VideoWriter('./videos/filter_panama.avi');
% writerObj.FrameRate = 5;
% open(writerObj);

figure('Renderer', 'painters', 'Position', [10 10 900 600])
for idx = 1:numel(warp.warp)
    idx
    frame_counter = frame_counter + round(numel(png_satellite_night.Files)/numel(warp.warp))-1;
    subtightplot(2,3,1)
    frameIndex = png_satellite_night.Files(frame_counter);
    satellite_image = imread(frameIndex{1,1});
    satellite_image_rot = imrotate(satellite_image,-90,'bilinear','crop');
    imshow(satellite_image_rot);axis off;title("Google Earth Night")
    subtightplot(2,3,2)
    N = histcounts2(warp.warp{1, idx}.y, warp.warp{1, idx}.x, pts, pts);
    [xG, yG] = meshgrid(-15:15);
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    imagesc(pts, pts, conv2(N, g, 'same'));axis equal;axis off;colormap("gray");title("Event-based Mosaicking")
%     imagesc(pts, pts, N);axis equal;axis off;colormap("gray");title("Event-based Mosaicking")
    subtightplot(2,3,[4 5])
    velocityX  = [velocityX;warp.warp{1, idx}.vx];
    velocityY  = [velocityY;warp.warp{1, idx}.vy];
    plot(velocityX,'r');hold on
    plot(velocityY,'b');
    legend("Vx","Vy");grid on;
    subtightplot(2,3,[3 6])
    imshow(panoramic_satellite);title("Panoramic Satellite Image")
    drawnow
    
%     pause(0.1)
%     F = getframe(gcf);
%     writeVideo(writerObj, F);
end
% close(writerObj);
% fprintf('Sucessfully generated the video\n')

% 2 plot one for warped event and one for heat map event side by side
% writerObj = VideoWriter('./videos/suez_canal_side_by_side_nofilter.avi');
% writerObj.FrameRate = 5;
% open(writerObj);
sigma = 0.5;
pts = linspace(-10, 200, 120);
for idx = 1:numel(warp.warp)
    figure(5675667)
    subtightplot(2,2,1)
    plot(abs(warp.warp{1, idx}.x), abs(warp.warp{1, idx}.y),'.');
    subtightplot(2,2,2)
    N = histcounts2(warp.warp{1, idx}.y, warp.warp{1, idx}.x, pts, pts);
    [xG, yG] = meshgrid(-15:15);
    g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
    g = g./sum(g(:));
    imagesc(pts, pts, conv2(N, g, 'same')*20);axis equal;axis off;colormap("gray");
    subplot(2,2,[3 4])
    velocityX  = [velocityX;warp.warp{1, idx}.vx];
    velocityY  = [velocityY;warp.warp{1, idx}.vy];
    plot(velocityX,'r');hold on
    plot(velocityY,'b');
    legend("Vx","Vy");grid on;
    drawnow
%     F = getframe(gcf);
%     writeVideo(writerObj, F);
end
% close(writerObj);
% fprintf('Sucessfully generated the video\n')


% subplot event warped patched
for idx = 1:numel(warp.warp)
    figure(5675667)
    subplot(9,9,idx)
    plot(abs(warp.warp{1, idx}.x), abs(warp.warp{1, idx}.y),'.');
end


% test
counter = 0;
time = 1e6:1e6:motion_compensated_image.ts(end);
for slice = time
    counter = counter + 1;
    tt = motion_compensated_image.ts;
    xx = motion_compensated_image.x;
    yy = motion_compensated_image.y;
    pp = motion_compensated_image.p;
    t0 = time(counter)+slice;
    iFirst = find(tt<t0,1,'last');
    iLast = find(tt>(t0+1e6),1,'first');
    ii = iFirst:iLast;
    disp((double(tt(iLast) - tt(iFirst)))/1e6)
    tt(ii);
    if ~isempty(ii)
        figure(5675667)
        subplot(9,9,counter)
        plot(motion_compensated_image.x(ii),motion_compensated_image.y(ii),'.');
    end
end
%% time surface
nEvents = numel(TD.x);
figure(2); clf;
S = zeros(max(TD.x),max(TD.y)); T = S; P = T;%-inf;
ss = imagesc(S); colorbar; axis image;
TD.ts = TD.ts(end) - TD.ts;
tau = 1e5;
displayFreq = 1; % in units of time
nextTimeSample = displayFreq;
for idx = 1:nEvents
    x = TD.x(idx);
    y = TD.y(idx);
    t = TD.ts(idx);
    p = TD.p(idx);
    T(x,y) = t;
    P(x,y) = p;
    if t > nextTimeSample
        nextTimeSample = nextTimeSample + displayFreq;
        S = P.*exp((T-t)/tau);
        set(ss,'CData',S)
        drawnow
    end
end
%% mosaicking
load('test_files/warp_filter_brittany.mat');
velocityX = [];
velocityY = [];
for idx = 1:numel(warp.warp)
    figure(567567)
    subplot(6,6,idx)
    plot(warp.warp{1, idx}.x, warp.warp{1, idx}.y,'.');
    velocityX  = [velocityX;warp.warp{1, idx}.vx];
    velocityY  = [velocityY;warp.warp{1, idx}.vy];
end
subplot(6,6,30)
sigma = 3;
pts = linspace(-2, 2, 50);
N = histcounts2(velocityX, velocityY, pts, pts);
[xG, yG] = meshgrid(-5:5);
g = exp(-xG.^2./(2.*sigma.^2)-yG.^2./(2.*sigma.^2));
g = g./sum(g(:));
imagesc(pts, pts, conv2(N, g, 'same')*20);xlabel("Vx");ylabel("Vy");
subplot(6,6,[31 36])
plot(velocityX,'r');hold on
plot(velocityY,'b');
legend("Vx","Vy");grid on;
%%
load('test_files/warp_filter_egypt.mat')
% load('/media/sam/Samsung_T5/PhD/Code/orbital_localisation/test_files/warp_nofilter_egypt.mat')

% combine all warp arrays into one single event stream
motion_compensated_imagex  = [];
motion_compensated_imagey  = [];
motion_compensated_imagep  = [];
motion_compensated_imagets = [];
motion_compensated_imagevx = [];
motion_compensated_imagevy = [];
for idx = 1:numel(warp.warp)
    motion_compensated_imagex  = [motion_compensated_imagex;warp.warp{1, idx}.x];
    motion_compensated_imagey  = [motion_compensated_imagey;warp.warp{1, idx}.y];
    motion_compensated_imagep  = [motion_compensated_imagep;warp.warp{1, idx}.p];
    motion_compensated_imagets = [motion_compensated_imagets;warp.warp{1, idx}.ts];
    motion_compensated_imagevx = [motion_compensated_imagevx;warp.warp{1, idx}.vx];
    motion_compensated_imagevy = [motion_compensated_imagevy;warp.warp{1, idx}.vx];
end
motion_compensated_image = struct("x",motion_compensated_imagex, ...
                                  "y",motion_compensated_imagey, ...
                                  "p",motion_compensated_imagep, ...
                                  "ts",motion_compensated_imagets, ...
                                  "vx",motion_compensated_imagevx, ...
                                  "vy",motion_compensated_imagevy);


xmin = min(motion_compensated_image.x);
ymin = min(motion_compensated_image.y);

xs = motion_compensated_image.x - xmin + 1;
ys = motion_compensated_image.y - ymin + 1;

pixels = zeros(700,700);

for idx = 1:numel(xs)
    pixels(ys(idx),xs(idx)) = (1.0 - xs(idx)) * (1.0 - ys(idx));
    pixels(ys(idx),xs(idx)+1) =  xs(idx) * (1.0 - ys(idx));
    pixels(ys(idx)+1,xs(idx)) =  1-xs(idx)*ys(idx);
    pixels(ys(idx)+1,xs(idx)+1) = xs(idx)*ys(idx);
end

figure(87878);
imagesc(pixels);colormap("gray");colorbar
