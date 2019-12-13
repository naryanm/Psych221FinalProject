%% Visual Acuity and Eye Movements Tool
% Final project script for PSYCH221, Fall 2019
% c. Naryan Murthy
% Based on script by Prof. Brian Wandell
% Requires ISETcam and ISETbio added to path to work. To use, modify
% values in User-Specified Parameters section, save, and hit 'Run.'

%% Clear and Initialize

clear all
close all
ieInit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% User-Specified Parameters

%%% Noise %%%
noiseType = 'frozen'; % Options: 'frozen'= no noise, 'random' = photon noise

%%% Eye Model %%%
modelType = 'control'; % Options:'center', 'periphery', 'control'
degEccentricity = 2; % For 'periphery' only, specify eccentricity in deg.

%%% Eye Movement Sequence %%%
nEyeMovements = 50; % integer # of movements
eyeMovementType = 'none'; 
%Options: 'none', 'horizontal', 'vertical', 'parallel', 'perpendicular', 'drift'

%For 'horizontal', 'vertical', and 'drift' only:
movementAmplitude = 8; %units = blocks of cone mosaic

%For 'horizontal' or 'vertical' only:
jitterAmplitude = 2; % units = blocks of cone mosaic

%%% Number of Trials %%%
nTrials = 1; %integer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Generate Slanted Edge Scene
%This never changes.

scene = sceneCreate('slanted edge');
oi = oiCreate;
oi = oiCompute(oi,scene);

%Uncomment to Display Scene 
%oiWindow(oi);
%sceneWindow(scene)
%% Initialize Cones

cones = coneMosaic;
cones.noiseFlag = noiseType;    %Set noise on or off
cones.integrationTime = 0.005;  %Five ms integration time

if strcmp(modelType, 'control')
    % Only L cones, standard cone aperture.
    cones.spatialDensity = [0,1,0,0]; 
elseif strcmp(modelType, 'center')
    % Only L and M cones, standard cone aperture.
    cones.spatialDensity = [0,0.5,0.5,0]; 
elseif strcmp(modelType, 'periphery')    
    %Set cone aperture to be bigger than standard, based on Curcio et al (1990).
    [spacing, aperture, density, params, comment] = coneSizeReadData('eccentricity',degEccentricity,'eccentricityUnits','deg','angle',degEccentricity);
    cones.pigment.pdHeight = aperture;
    cones.pigment.pdWidth = aperture;    
    %Use L, M, and S cones according to ratio in Foundations of Vision
    cones.spatialDensity = [0,0.6,0.3,0.1];
else   
    error('Unrecognized model type.')    
end


% Loop through all trials, computing and storing MTF.
for j=1:nTrials
%% Create Eye Movements (Horizontal or Vertical)
    
MovementPositions = zeros(nEyeMovements,2);
[x,y]=size(cones.pattern);

if strcmp(eyeMovementType, 'vertical')
    for i=1:nEyeMovements
        if isodd(i)
            MovementPositions(i,1)=jitterAmplitude*round(round(1-2*rand(1)));
            MovementPositions(i,2)=-movementAmplitude;
        else
            MovementPositions(i,1)=jitterAmplitude*round(round(1-2*rand(1)));
            MovementPositions(i,2)=movementAmplitude;
        end
    end
    %First movement needs to be [0,0] to locate edge of scene later.
    MovementPositions(1,:)=[0,0];
    cones.emPositions = MovementPositions;
elseif strcmp(eyeMovementType, 'horizontal')
    for i=1:nEyeMovements
        if isodd(i)
            MovementPositions(i,2)=jitterAmplitude*round(round(1-2*rand(1)));
            MovementPositions(i,1)=-movementAmplitude;
        else
            MovementPositions(i,2)=jitterAmplitude*round(round(1-2*rand(1)));
            MovementPositions(i,1)=movementAmplitude;
        end
    end
    %First movement needs to be [0,0] to locate edge of scene later.
    MovementPositions(1,:)=[0,0];
    cones.emPositions = MovementPositions;
elseif strcmp(eyeMovementType, 'parallel')
    for i=1:nEyeMovements
        if isodd(i)
            MovementPositions(i,1)=round(x/48);
            MovementPositions(i,2)=-round(y/24);
        else
            MovementPositions(i,1)=-round(x/48);
            MovementPositions(i,2)=round(y/24);
        end
    end
    %First movement needs to be [0,0] to locate edge of scene later.
    MovementPositions(1,:)=[0,0];
    cones.emPositions = MovementPositions;
elseif strcmp(eyeMovementType, 'perpendicular')
    for i=1:nEyeMovements
        if isodd(i)
            MovementPositions(i,1)=-round(y/24);
            MovementPositions(i,2)=-round(x/48);
        else
            MovementPositions(i,1)=round(y/24);
            MovementPositions(i,2)=round(x/48);
        end
    end
    %First movement needs to be [0,0] to locate edge of scene later.
    MovementPositions(1,:)=[0,0];
    cones.emPositions = MovementPositions;
elseif strcmp(eyeMovementType, 'fixational')
    em = fixationalEM();
    cones.emGenSequence(nEyeMovements,'rSeed',[]); % Randomize eye movements
    cones.compute(oi);
elseif strcmp(eyeMovementType, 'drift')
    
    MovementPositions(1:floor(nEyeMovements/4),1) = linspace(0,-movementAmplitude,floor(nEyeMovements/4));
    
    MovementPositions((floor(nEyeMovements/4)+1):2*floor(nEyeMovements/4),1) = linspace(-movementAmplitude, 0, floor(nEyeMovements/4));
           
    MovementPositions((2*floor(nEyeMovements/4)+1):3*floor(nEyeMovements/4),1) = linspace(0,movementAmplitude,floor(nEyeMovements/4));
  
    MovementPositions((3*floor(nEyeMovements/4)+1):4*floor(nEyeMovements/4),1) = linspace(movementAmplitude, 0,floor(nEyeMovements/4));

    %First movement needs to be [0,0] to locate edge of scene later.
    MovementPositions(1,:)=[0,0];
    cones.emPositions = MovementPositions;
    
elseif strcmp(eyeMovementType, 'none')
    cones.emPositions = zeros(nEyeMovements,2);
else
    error('Unrecognized Eye Movement Type')
end

%% Visualize Eye Movement Pattern

%Display all fixational eye movements, mostly out of curiosity. 
%This only runs if the number of trials is less than or equal to 10, to keep the number
%of plots generated reasonable. 

if (nTrials <= 10) && strcmp(eyeMovementType,'fixational')
figure
plot(cones.emPositions(:,1), cones.emPositions(:,2),'r.-')
axis([-88/2, 88/2, -72/2, 72/2])
title([eyeMovementType,' Eye Movement Pattern, Trial ',num2str(j)])
xlabel('Horizontal Position (cones)')
ylabel('Vertical Position (cones)')
end

%% Compute and Store MTF Info

    cones.compute(oi);
    absorptions = cones.absorptions;
    edgeImage  = absorptions(:,:,1); %True Edge Image is from first frame.

    % ISO12233 requires cropping the image so that the edge is taller than wide
    % This is a specific way to find the rect for the image
    rect = ISOFindSlantedBar(edgeImage);
      %ieNewGraphWin; imagesc(edgeImage); colormap(gray); axis image
      %h = drawrectangle('Position',rect);

    % Now the whole set of eye movements
     edgeImage = mean(absorptions,3);
     edgeImage = imcrop(edgeImage,rect);
     %ieNewGraphWin; imagesc(edgeImage); axis image; colormap(gray);

    % Call the ISO routine
    % It returns the key parameters
    deltaX = cones.patternSampleSize(1)*1e3;
    mtf{j} = ISO12233(edgeImage,deltaX,[],'none');
    
    %Truncate below Nyquist freq.
    mtf{j}.mtf = mtf{j}.mtf(mtf{j}.freq < mtf{j}.nyquistf);
    mtf{j}.freq = mtf{j}.freq(mtf{j}.freq < mtf{j}.nyquistf);

end

%% Visualize Cone Representation of Image & Cone Mosaic
%This is for last trial only
cones.window;

%% Plot MTF50 of All Trials

MTF50 = NaN*ones(1,nTrials);
mmPerDeg = 0.2852;

for k = 1:nTrials
    if isfield(mtf{k},'mtf50') 
        %Only keep MTF50 value if it's below Nyquist freq. and is less than
        %or equal to 15 cycles/degree.
        if (mtf{k}.mtf50 < mtf{k}.nyquistf) && (mtf{k}.mtf50*mmPerDeg <= 15) && (mtf{k}.mtf50*mmPerDeg > 0)
            MTF50(k) = mtf{k}.mtf50*mmPerDeg;
        else
            MTF50(k) = NaN;
        end
    end
end

figure
plot(1:nTrials, MTF50, 'b*-');
axis([0 nTrials 0 20])
title(['MTF50 for All ',num2str(nTrials),' Trials, Excluding Erroneous Points'])
xlabel('Trial Number')
ylabel('MTF50 Spatial Frequency (cycles/degree)')

%Compute Mean of valid data points, and print to command window.
MTF50_mean = mean(MTF50(~isnan(MTF50)))
MTF50_stdev = std(MTF50(~isnan(MTF50)))

%% Plot MTFs of All Trials

figure
for m=1:nTrials
    plot(mtf{m}.freq*mmPerDeg, mtf{m}.mtf)
    xlabel('Spatial Frequency (cycles/deg)');
ylabel('Contrast Reduction (SFR)');
grid on;
axis([0 60 0 1])
title(['All MTFs for ',modelType,' Eye Model, ',num2str(nEyeMovements),...
        ' ',eyeMovementType, ' Eye Movements, and ',num2str(nTrials),' Trial(s)'])
    hold on
end
