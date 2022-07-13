%% Generate dataset (Fake star moving across the sky)
number_of_stars = 1;
ys = 500;
xs = ys;
R = 20;
speed = 1;
noise = 0.001;
I0 = zeros(ys,xs);

e = []; idx = 0;
for star = 1:number_of_stars
    randomx = randperm(size(I0,1),1);
    randomy = randperm(size(I0,1),1);
    centre =  [randomx,randomy]; %[1+R,1+R];
%     centre = [25 1];
    [X,Y]=ndgrid(1:size(I0,1),1:size(I0,2));
    X=X-centre(1);Y=Y-centre(2);
    L=sqrt(X.^2+Y.^2)<=R;
    
    for shift = 1:speed:ys-(R*2+2)
        II = imtranslate(L,[shift, 1]);
        J = imnoise(double(II),'salt & pepper',noise);
        for y=1:ys
            for x = 1:xs
                if J(x,y)~= 0 && rand >((ys-1)/ys)
                    idx = idx + 1;
                    e.x(idx) = x;
                    e.y(idx) = y;
                    e.p(idx) = J(x,y);
                    e.t(idx) = idx;
                end
            end
        end
    end
end
figure(6786868);
scatter3(e.x,e.y,e.t,'.');

% nEvents = numel(e.x);
% figure(2); clf;
% S = zeros(xs,ys); T = S; P = T;%-inf;
% ss = imagesc(S); colorbar; axis image;

% tau = 2e2;
% displayFreq = 1; % in units of time
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

counter = 0;
VelX = [];
VelY = [];
warped_events_x = [];
warped_events_y = [];

range = 2e2;
timeArray = range:range:e.t(end);
for shunck = timeArray
    counter = counter + 1;
    counter
    ii = find(e.t>shunck & e.t<(shunck+range));
    
    if ~isempty(ii)
        for iVelx = 1:nVel
            vx = VelArray(iVelx);
            for x = 1:xs
                xSum(x) = sum(round(e.x(ii)+vx*e.t(ii)) == x);
            end
            xStd(iVelx) = std(xSum);
        end
        
        for iVely = 1:nVel
            vy = VelArray(iVely);
            for y = 1:ys
                ySum(y) = sum(round(e.y(ii)+vy*e.t(ii)) == y);
            end
            yStd(iVely) = std(ySum);
        end
        
        [Vx,indx] = max(xStd);
        [Vy,indy] = max(yStd);
        
        velocity_x = VelArray(indx);
        velocity_y = VelArray(indy);
        VelX = [VelX;velocity_x];
        VelY = [VelY;velocity_y];
        warped_events_x = [warped_events_x;round(e.x(ii)+velocity_x*e.t(ii))'];
        warped_events_y = [warped_events_y;round(e.y(ii)+velocity_y*e.t(ii))'];
    end
end
e.warpedx  = warped_events_x;
e.warpedy  = warped_events_y;

figure(68678);
subplot(2,4,1)
scatter3(e.x,e.y,e.t,'.');title("Input Events");
subplot(2,4,2)
Y_x = (1:numel(VelX))';
Y_y = (1:numel(VelY))';
plot_shaded(Y_x,normalize(VelX,'range'));hold on;plot_shaded(Y_y,normalize(VelY,'range'));
title("Optical Velocity \color{orange} Vx, \color{blue} Vy");grid on;
subplot(2,4,[3 4])
scatter3(e.warpedx,e.warpedy,1:numel(e.warpedy),'.');title("3D Warped Events");
subplot(2,4,[5 6])
plot(e.warpedx,e.warpedy,'.');title("2D Warped Events");
subplot(2,4,[7 8])
sigma = 0.7;
pts = linspace(min(e.warpedx), max(e.warpedx), 80);
N = histcounts2(e.warpedx, e.warpedy, pts, pts);
imagesc(pts, pts, N);colorbar;title("Warped Events HeatMap");
