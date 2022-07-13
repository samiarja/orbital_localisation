%% Generate dataset (Fake star moving across the sky)
number_of_stars = 2;
ys = 500;
xs = ys;
R = 5;
speed = 1;
noise = 0.001;
I0 = zeros(ys,xs);

e = []; idx = 0;
for star = 1:number_of_stars
    randomx = randperm(size(I0,1),1);
    randomy = randperm(size(I0,1),1);
    centre = [randomx,randomy]; %[1+R,1+R];
    [X,Y]=ndgrid(1:size(I0,1),1:size(I0,2));
    X=X-centre(1);Y=Y-centre(2);
    L=sqrt(X.^2+Y.^2)<=R;
    for shift = 1:speed:ys-(R*2+2)
        II = imtranslate(L,[shift, shift]);
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
% 
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
VelArray = [-10:.01:10];
nVel = numel(VelArray);
xStd = nan(1,nVel);
xSum = nan(1,xs);
yStd = nan(1,nVel);
ySum = nan(1,ys);

for iVelx = 1:nVel
    vx = VelArray(iVelx);
    for x = 1:xs
        xSum(x) = sum(round(e.x+vx*e.t) == x);
    end
    xStd(iVelx) = std(xSum);
end

for iVely = 1:nVel
    vy = VelArray(iVely);
    for y = 1:ys
        ySum(y) = sum(round(e.y+vy*e.t) == y);
    end
    yStd(iVely) = std(ySum);
end

[Vx,indx] = max(xStd);
[Vy,indy] = max(yStd);

velocity_x = VelArray(indx);
velocity_y = VelArray(indy);

e.warpedx  = round(e.x+velocity_x*e.t);
e.warpedy  = round(e.y+velocity_y*e.t);

figure(68678);
subplot(2,3,1)
scatter3(e.x,e.y,e.t,'.');title("Input Events");
subplot(2,3,2)
plot(xStd,yStd,'.');hold on
plot(Vx,Vy,'or');grid on;title("Std for X and Y");
subplot(2,3,3);
plot(xStd,'r');hold on
plot(yStd,'b');grid on;
plot(indx,Vx,'or');
plot(indy,Vy,'ob');title("Vx and Vy");
legend("xStd","yStd","Vx","Vy");
subplot(2,3,4)
scatter3(e.warpedx,e.warpedy,e.t,'.');title("3D Warped Events");
subplot(2,3,5)
plot(e.warpedx,e.warpedy,'.');title("2D Warped Events");
subplot(2,3,6)
sigma = 0.7;
pts = linspace(min(e.warpedx), max(e.warpedx), 80);
N = histcounts2(e.warpedx, e.warpedy, pts, pts);
imagesc(pts, pts, N);colorbar;title("Warped Events HeatMap");
