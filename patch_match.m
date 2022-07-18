%% stitch by patch single frame
event     = imread('data/Egypt/egypt_events_cropped_suez_canal.png');
satellite = imread('data/Egypt/egypt_satellite_cropped_suez_canal.png');

eventImage = rgb2gray(satellite);
mapImage   = rgb2gray(event);

% detect image feature points
boxPoints   = detectSURFFeatures(eventImage);
scenePoints = detectSURFFeatures(mapImage);

figure(676767);
subtightplot(1,4,1)
imshow(eventImage);
title('Event');
hold on;
plot(selectStrongest(boxPoints, 10));

subtightplot(1,4,2)
imshow(mapImage);
title('Map');
hold on;
plot(selectStrongest(scenePoints, 10));

% Extract feature descriptor
[boxFeatures, boxPoints] = extractFeatures(eventImage, boxPoints);
[sceneFeatures, scenePoints] = extractFeatures(mapImage, scenePoints);

% find Putative point matches
boxPairs = matchFeatures(boxFeatures, sceneFeatures);

% display matched features
matchedBoxPoints = boxPoints(boxPairs(:, 1), :);
matchedScenePoints = scenePoints(boxPairs(:, 2), :);
subtightplot(1,4,3)
showMatchedFeatures(eventImage, mapImage, matchedBoxPoints, matchedScenePoints, 'montage');
title('Putatively Matched Points(Including Outliers)');

% Locate the Object in Scene Using Putative Matches
[tform, inlierBoxPoints, inlierScenePoints] = ...
    estimateGeometricTransform(matchedBoxPoints, matchedScenePoints, 'affine');
% subtightplot(1,4,4);
% showMatchedFeatures(eventImage, mapImage, inlierBoxPoints, inlierScenePoints, 'montage');
% title('Matched Points(Inliers Only)');

% Get the bounding polygon of the reference image
boxPolygon = [1, 1;...                              %top-left
        size(eventImage, 2), 1;...                    %top-right
        size(eventImage, 2), size(eventImage, 1);...    %bottom-right
        1, size(eventImage, 1);...                    %bottom-left
        1, 1];                                      %top-left again to close the polygon

% Transform the polygon into the coordinate system of the target image
% The transformed polygon indicates the location of the object in scene
newBoxPolygon = transformPointsForward(tform, boxPolygon);
subtightplot(1,4,4)
imshow(mapImage);
hold on;
line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'y');
title('Detected Box');

%% stitch by patch video
event     = imageDatastore('data/Egypt/png_event/');
satellite = imageDatastore('data/Egypt/png_satellite/');

for idx = 1:numel(satellite.Files)
    satelliteFiles = satellite.Files(idx);
    boxImage = imread(satelliteFiles{1,1});
    boxImage = rgb2gray(boxImage);
    
    satelliteFiles = event.Files(idx);
    sceneImage = imread(satelliteFiles{1,1});
    sceneImage = rgb2gray(sceneImage);
    
    % detect image feature points
    boxPoints = detectSURFFeatures(boxImage);
    scenePoints = detectSURFFeatures(sceneImage);
    figure(676767);
    subtightplot(2,2,1)
    imshow(boxImage);
    title('100 Feature Points');
    hold on;
    plot(selectStrongest(boxPoints, 100));
    
    subtightplot(2,2,2)
    imshow(sceneImage);
    title('300 Feature Points');
    hold on;
    plot(selectStrongest(scenePoints, 300));
    
    % Extract feature descriptor
    [boxFeatures, boxPoints] = extractFeatures(boxImage, boxPoints);
    [sceneFeatures, scenePoints] = extractFeatures(sceneImage, scenePoints);
    
    % find Putative point matches
    boxPairs = matchFeatures(boxFeatures, sceneFeatures);
    
    % display matched features
    matchedBoxPoints = boxPoints(boxPairs(:, 1), :);
    matchedScenePoints = scenePoints(boxPairs(:, 2), :);
    subtightplot(2,2,[3 4])
    showMatchedFeatures(boxImage, sceneImage, matchedBoxPoints, matchedScenePoints, 'montage');
    title('Putatively Matched Points(Including Outliers)');
    
end

% % % % % % Locate the Object in Scene Using Putative Matches
% % % % % [tform, inlierBoxPoints, inlierScenePoints] = ...
% % % % %     estimateGeometricTransform(matchedBoxPoints, matchedScenePoints, 'affine');
% % % % % figure(4);
% % % % % showMatchedFeatures(boxImage, sceneImage, inlierBoxPoints, inlierScenePoints, 'montage');
% % % % % title('Matched Points(Inliers Only)');
% % % % %
% % % % % % Get the bounding polygon of the reference image
% % % % % boxPolygon = [1, 1;...                              %top-left
% % % % %         size(boxImage, 2), 1;...                    %top-right
% % % % %         size(boxImage, 2), size(boxImage, 1);...    %bottom-right
% % % % %         1, size(boxImage, 1);...                    %bottom-left
% % % % %         1, 1];                                      %top-left again to close the polygon
% % % % %
% % % % % % Transform the polygon into the coordinate system of the target image
% % % % % % The transformed polygon indicates the location of the object in scene
% % % % % newBoxPolygon = transformPointsForward(tform, boxPolygon);
% % % % % figure(5);
% % % % % imshow(sceneImage);
% % % % % hold on;
% % % % % line(newBoxPolygon(:, 1), newBoxPolygon(:, 2), 'Color', 'y');
% % % % % title('Detected Box');

%% Registering to a Digital Orthophoto (Egypt Slice 1)
inputMapEventsUnregistered = imread('data/Egypt/egypt_event_slice1.png');
baseSatelliteImage         = imread('data/Egypt/egypt_satellite_slice1.png');
load('data/Egypt/cpstruct_slice1.mat')

% cpselect(baseSatelliteImage(:,:,1),inputMapEventsUnregistered);

mytform = cp2tform(cpstruct.inputPoints,cpstruct.basePoints,'projective');
baseSatelliteImage_transformed = imtransform(baseSatelliteImage,mytform);
figure(678678);
II = baseSatelliteImage_transformed;
J = imtranslate(II,[10, -14]);
subtightplot(1,4,1)
image(J);title("Registered Satellite map of Suez Canal (Night)");axis off
subtightplot(1,4,2)
image(inputMapEventsUnregistered);title("Image of the Warped Event (IWP)");axis off
subtightplot(1,4,3)
image(J);
axis off
hold on
im = image(inputMapEventsUnregistered);
im.AlphaData = max(inputMapEventsUnregistered,[],3)+30;axis off
hold off;
title("Overlap of the satellite map and the IWP");
subtightplot(1,4,4)

Jnew = imtranslate(J,[20, 15]);
paddedJnew = padarray(Jnew,[31 0]);paddedJnew(end,:,:) = [];
paddedinputMap = padarray(inputMapEventsUnregistered,[4 10]);paddedinputMap(:,end,:) = [];

[ssimval,ssimmap] = ssim(paddedinputMap,paddedJnew);
imshow(ssimmap,[]);title("Structural similarity (SSIM)");axis normal

%% Registering to a Digital Orthophoto (Egypt Slice 2)
inputMapEventsUnregistered = imread('data/Egypt/egypt_event_slice2.png');
baseSatelliteImage         = imread('data/Egypt/egypt_satellite_slice2.png');
% load('data/Egypt/cpstruct_slice1.mat')

cpselect(baseSatelliteImage(:,:,1),inputMapEventsUnregistered);

mytform = cp2tform(cpstruct3.inputPoints,cpstruct3.basePoints,'projective');
baseSatelliteImage_transformed = imtransform(baseSatelliteImage,mytform);
figure(678678);
II = baseSatelliteImage_transformed;
J = imtranslate(II,[10, -14]);
subtightplot(1,4,1)
image(J);title("Registered Satellite map of Suez Canal (Night)");axis off
subtightplot(1,4,2)
image(inputMapEventsUnregistered);title("Image of the Warped Event (IWP)");axis off
subtightplot(1,4,3)
image(J);
axis off
hold on
im = image(inputMapEventsUnregistered);
im.AlphaData = max(inputMapEventsUnregistered,[],3)+30;axis off
hold off;
title("Overlap of the satellite map and the IWP");
subtightplot(1,4,4)

Jnew = imtranslate(J,[20, 15]);
paddedJnew = padarray(Jnew,[31 0]);paddedJnew(end,:,:) = [];
paddedinputMap = padarray(inputMapEventsUnregistered,[4 10]);paddedinputMap(:,end,:) = [];

[ssimval,ssimmap] = ssim(paddedinputMap,paddedJnew);
imshow(ssimmap,[]);title("Structural similarity (SSIM)");axis normal

%% Registering to a Digital Orthophoto (Egypt Slice 3)
inputMapEventsUnregistered = imread('data/Egypt/egypt_event_slice3.png');
baseSatelliteImage         = imread('data/Egypt/Screenshot from 2022-07-18 11-45-35.png');
% load('data/Egypt/cpstruct_slice1.mat')

% cpselect(baseSatelliteImage(:,:,1),inputMapEventsUnregistered);
% 
% mytform = cp2tform(cpstruct1.inputPoints,cpstruct1.basePoints,'projective');
% baseSatelliteImage_transformed = imtransform(baseSatelliteImage,mytform);
% figure(678678);
% II = baseSatelliteImage_transformed;
% J = imtranslate(II,[10, -14]);
% subtightplot(1,4,1)
% image(J);title("Registered Satellite map of Suez Canal (Night)");axis off
% subtightplot(1,4,2)
% image(inputMapEventsUnregistered);title("Image of the Warped Event (IWP)");axis off
% subtightplot(1,4,3)
% image(J);
% axis off
% hold on
% im = image(inputMapEventsUnregistered);
% im.AlphaData = max(inputMapEventsUnregistered,[],3)+30;axis off
% hold off;
% title("Overlap of the satellite map and the IWP");
% subtightplot(1,4,4)
% 
% Jnew = imtranslate(J,[20, 15]);
% paddedJnew = padarray(Jnew,[31 0]);paddedJnew(end,:,:) = [];
% paddedinputMap = padarray(inputMapEventsUnregistered,[4 10]);paddedinputMap(:,end,:) = [];
% 
% [ssimval,ssimmap] = ssim(paddedinputMap,paddedJnew);
% imshow(ssimmap,[]);title("Structural similarity (SSIM)");axis normal


inputMapEventsUnregistered = imread('data/Egypt/egypt_event_slice3.png');
baseSatelliteImage         = imread('data/Egypt/Screenshot from 2022-07-18 11-45-35.png');
figure(678686);
subtightplot(1,4,3)
eventJ = imrotate(inputMapEventsUnregistered,5);
image(eventJ);
axis off
hold on
satelliteJ = imtranslate(baseSatelliteImage,[-35, 5]);
im = image(satelliteJ);
im.AlphaData = max(baseSatelliteImage,[],3)+70;axis off
hold off;axis normal;title("Overlap of the satellite map and the IWP");
subtightplot(1,4,1)
image(satelliteJ);axis off;title("Registered Satellite map of Suez Canal (Night)");
subtightplot(1,4,2)
image(eventJ);
image(eventJ);axis off;title("Image of the Warped Event (IWP)");
subtightplot(1,4,4)
eventJnew = padarray(eventJ,[0 18]);eventJnew(:,end,:) = [];
satelliteJnew = padarray(satelliteJ,[53 1]);satelliteJnew(:,end,:) = [];
[ssimval,ssimmap] = ssim(eventJnew,satelliteJnew);
imshow(ssimmap,[]);title("Structural similarity (SSIM)");axis normal

%% Brittany registration
inputMapEventsUnregistered = imread('data/brittany/event.png');
baseSatelliteImage         = imread('data/brittany/satellite2.png');

inputMapEventsUnregisteredRot   = imrotate(inputMapEventsUnregistered,270);
baseSatelliteImageRot           = imrotate(baseSatelliteImage,-6);
% figure(6778);
% montage({inputMapEventsUnregisteredRot,baseSatelliteImageRot});
% 
% cpselect(baseSatelliteImageRot(:,:,1),inputMapEventsUnregisteredRot);
% mytform = cp2tform(cpstruct1.inputPoints,cpstruct1.basePoints,'projective');
% baseSatelliteImage_transformed = imtransform(baseSatelliteImageRot,mytform);

figure(678678);
% J = imtranslate(II,[10, -14]);
subtightplot(1,4,1)
% baseSatelliteImageRot = imrotate(baseSatelliteImage,260);
image(baseSatelliteImageRot);title("Registered Satellite map of Brittany (Night)");axis off
subtightplot(1,4,2)
image(inputMapEventsUnregisteredRot);title("Image of the Warped Event (IWP)");axis off
subtightplot(1,4,3)
image(baseSatelliteImageRot);axis off;hold on
inputMapEventsUnregisteredRotResized = imresize(inputMapEventsUnregisteredRot,0.75);
eventShifted = imtranslate(inputMapEventsUnregisteredRotResized,[-10, -15]);
eventShiftedRot = imrotate(eventShifted,-10);
im = image(eventShiftedRot);
im.AlphaData = max(eventShiftedRot,[],3)+30;axis off
hold off;
title("Overlap of the satellite map and the IWP");

subtightplot(1,4,4)
paddedJnew = padarray(baseSatelliteImageRot,[8 0]);paddedJnew(end,:,:) = [];
paddedinputMap = padarray(eventShiftedRot,[0 22]);%paddedinputMap(:,end,:) = [];

[ssimval,ssimmap] = ssim(paddedinputMap,paddedJnew);
imshow(ssimmap,[]);title("Structural similarity (SSIM)");axis normal

%% Panama registration
inputMapEventsUnregistered = imread('data/panama/Panama_2022-01-24_20_12_11_NADIR_21.49_-0.74.png');
baseSatelliteImage         = imread('data/panama/satellite.png');
inputMapEventsUnregisteredRot           = imrotate(inputMapEventsUnregistered,272);

% figure(456567)
% montage({baseSatelliteImage,inputMapEventsUnregisteredRot})

% cpselect(baseSatelliteImage(:,:,1),inputMapEventsUnregisteredRot);
% mytform = cp2tform(cpstruct.inputPoints,cpstruct.basePoints,'projective');
% baseSatelliteImage_transformed = imtransform(baseSatelliteImage,mytform);

figure(678678);
II = baseSatelliteImage_transformed;
J = imtranslate(II,[10, -14]);
subtightplot(1,4,1)
image(baseSatelliteImage);
title("Registered Satellite map of Panama (Day)");axis off
subtightplot(1,4,2)
image(inputMapEventsUnregisteredRot);
title("Image of the Warped Event (IWP)");axis off
subtightplot(1,4,3)
image(baseSatelliteImage);axis off;hold on
eventShifted = imtranslate(inputMapEventsUnregisteredRot,[20, -52]);
im = image(eventShifted);
im.AlphaData = max(eventShifted,[],3)+30;axis off
hold off;
title("Overlap of the satellite map and the IWP");
subtightplot(1,4,4)
% paddedEvent     = imtranslate(eventShifted,[20, 15]);
paddedSatellite = padarray(baseSatelliteImage,[55 5]);paddedSatellite(:,end,:) = [];
paddedinputMap  = padarray(inputMapEventsUnregisteredRot,[4 10]);paddedinputMap(:,end,:) = [];

[ssimval,ssimmap] = ssim(paddedSatellite,eventShifted);
imshow(ssimmap,[]);title("Structural similarity (SSIM)");axis normal