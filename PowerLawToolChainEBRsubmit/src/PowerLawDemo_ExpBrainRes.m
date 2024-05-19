% Demo Code for One-Third Power Law Principled Protocol
%
% Run this script to generate synthetic ellipse trajectories
% obeying the 1/3 power law, or linear velocites as per
% Maoz et al., 2005 DOI: 10.13140/2.1.4401.8884
% Schaal & Sternad, 2001 DOI: 10.1007/s002210000505
%
% white / pink noise of a range of Standard Deviations may be added
%
% Velocity Gain Factor and Beta Power Law exponents are then calculated via
% choices of differentiation, filtering, curvature calculation, and regression
% as outlined in Exp Brain Res review paper of Fraser et al., 2024
%
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
% requires functions subfolder
% - curvatureKinematicEBR.m
% - differentiateKinematicsEBR.m
% - curvatureMengerEBR.m
% - regressDataEBR.m
%
% require req subfolder with
% https://uk.mathworks.com/matlabcentral/fileexchange/34874-interparc
% - interparc/
%             -interparc.m
%

close all
clear all
commandwindow

% add the subfunction in these folders
addpath(genpath('functions'))
addpath(genpath('req'))

%% control display of intermediary and saving of final images.
saveImage = 1; % save the final images of Velocity Gain Function and Beta Exponent against Noise
displayFilters = 0; % display intermediate output of diffs + SG & LPBW
displayGraphs = 0; % display intermediate output for regressions
displayCurvature = 0; % display delta between Menger and Kinematic Curvature
% intermediate graphs are not saved

%% analysis choices
% linear velocity or 1/3 power law compliant ellipses
equidistant = 1; % 2 for Linear tangential velocity ellipse trajectory
% 1 for 1/3 power law tangential velocity

% paramChoice - options to generate synthetic data matching
paramChoice = 1; % 1 for Maoz, 2 for Schaal
Maoz = 1; % Maoz et al., 2005 DOI: 10.13140/2.1.4401.8884
Schaal = 2; % Schaal & Sternad, 2001 DOI: 10.1007/s002210000505
expStr = {'Maoz' 'Schaal'};
% if the data is 1/3 power law compliant, then only the geometry is used

% curvature calculation choice
% geometric Menger curvature versus Curvature from Velocity and Acceleration
curvatureChoice = 1; % 1 for Kinematic Curvature, 2 for Menger Curvature
kinematicCurvature = 1; %
mengerCurvature = 2; %
kinStr = {'KinematicCurvature' 'GeometricCurvature'};

edgeClip = 50; % differentiation and Savitzky-Golay have edge effects
% we clip the edges when generating the kinematics from the raw displacement.
% edgeClip should be greater than the width of the SG filter.

% strings for graphs and filenames
noiseStr = {'White' 'Pink'};
filterStr = {'Diff' 'Diff-LPBW' 'LPBW-Diff' 'S-G'};
regressStr = {'regress0' 'fitlm0' 'fitlmY' 'LMLS' 'IRLS' };
equiStr = {'OneThirdPowerLaw' 'LinearMotion'};

% noise generation
stepSize = 1; % noise Standard Deviations, 0 to maxVal with steps between
maxVal = 5; % max noise Std Dev in mm
reps = 10; % number of unique noise generations at each Std Dev, higher gives smoother graph
noiseType = 1; % 1 is white noise, 2 is pink (which is ~30 times smaller in Std when generated so we multiple up)
pinkNoiseMod = 30; % pink noise is ~30x less Std Deviation, this rectifies.
noise = [];
noiseLength = length([0:stepSize:maxVal]);
% generate noise at these Standard Deviations
for rep = 1:reps
    noise = [noise [0:stepSize:maxVal]];
end

%% Generate ellipses
switch paramChoice
    case 1 % Maoz
        disp('generating Maoz et al., 2005, style ellipses')
        filterParamsBW = [2 10]; % Butterworth      [ Order - Low Pass Corner Freq ]
        filterParamsSG = [4 17]; % Savitzky-Golay   [ Order - Width ]
        % As per Crenna et al., 2021 DOI 1424-8220/21/13/4580
        numPoints = 120; % number of linearly spaced points around the ellipse to generate
        Fc = 0.5; % 1/period = frequency of ellipse traces per second
        majAxis = 350; % ellipse major axis in mm
        minAxis = 130; % ellipse minor axis in mm
    case 2 %Schaal
        disp('generating Schaal and Sternad, 2001, style ellipses')
        filterParamsBW = [4 2.5]; % [ Order - Low Pass Corner Freq ]
        filterParamsSG = [4 17]; % Savitzky-Golay   [ Order - Width ]
        % As per Crenna et al., 2021 DOI 1424-8220/21/13/4580
        numPoints = 6;
        Fc = 0.5;
        majAxis = 200;
        minAxis = 100;
end

%% Beta and VGF (yGain) parameter guesses for nonlinear regressions
betaSeed = -1/3;
yGainSeed = 0.5;
LMSeeds = [yGainSeed betaSeed];

% loop over the length of the Std Deviation noise array
% generating an ellipse (x 'reps' times) and then running
% legacy and principled analysis protocols

for noises = 1:length(noise)

    dt = 0.01;
    fs = 1/dt;
    t = (0:dt:10)';
    noiseStd = noise(noises); % the current noise for this loop

    % compose an Ellipse with sin in the x-axis and cos in the Y
    % contaminated with noise

    xPoints = (majAxis/2)*sin(2*pi*Fc*t); % mm main axis ellipse
    yPoints = (minAxis/2)*cos(2*pi*Fc*t); % mm minor axis ellipse
    % this composes ellipses but the points are not equidistant.
    % we use https://uk.mathworks.com/matlabcentral/fileexchange/34874-interparc
    % to get equidistant points 

    if equidistant == 2 % we take linear velocity points along the hull

        pt = interparc(numPoints*5,xPoints,yPoints,'spline');
        % get points along the hull of each of the generated ellipse equidistantly
        % when the ellipse Fc is 0.5, we get 5 full ellipse traces in 10
        % seconds

        xPure = pt(:,1);
        yPure = pt(:,2);

        % generate noise for both x and y of suitable length, in white or
        % pink
        if noiseType == 1 % white

            xNoise = noiseStd.*randn(size(xPure)); % add noise with Standard Deviation
            yNoise = noiseStd.*randn(size(yPure));

        elseif noiseType == 2 % pink

            xNoise = pinkNoiseMod.*noiseStd.*pinknoise(size(xPure)); % add noise with Standard Deviation
            yNoise = pinkNoiseMod.*noiseStd.*pinknoise(size(yPure));

        end

    else % we keep the sin and cos, which perfectly reproduce the one-third power law

        xPure = xPoints;
        yPure = yPoints;

        if noiseType == 1 % white

            xNoise = noiseStd.*randn(size(xPoints)); % add noise with Standard Deviation
            yNoise = noiseStd.*randn(size(yPoints));

        elseif noiseType == 2 % pink

            xNoise = pinkNoiseMod.*noiseStd.*pinknoise(size(xPoints)); % add noise with Standard Deviation
            yNoise = pinkNoiseMod.*noiseStd.*pinknoise(size(yPoints));

        end

    end

    figNum= 9;

    % we combine the linear OR one-third power law compliant ellipse
    % trajectory with the white or pink noise of a specific Std Deviation.
    % given the displacements in x and y, we generate velocity and
    % acceleration to give curvature

    kinematicDataX = [];
    kinematicDataY = [];

    for noiseTest = 2
        if noiseTest == 2
            x = xPure + xNoise;
            y = yPure + yNoise;
            noisy  = 'Noisy';
        else
            x = xPure;
            y = yPure;
            noisy = 'Pure';
        end
        for filterType = 1:4

            if filterType == 4 % SG
                [dx, dy]= differentiateKinematicsEBR(x, y, filterType, filterParamsSG, fs);
            else
                [dx, dy]= differentiateKinematicsEBR(x, y, filterType, filterParamsBW, fs);
            end

            drawnow
            figNum = figNum +1;
            if displayFilters
                figure(figNum)

                plot(x(20:end-10),'.-')
                hold on
                plot(dx(20:end-10,2)) %velocity
                plot(dx(20:end-10,3)) %acceleration
                %plot(dx(20:end-10,4)) % jerk
                hold off

                %         xlim([0 900 ]);
                %         ylim([-4000 4000]);

                %% calculate beta - we now have acceleration and velocity!

                legend('x','x`','x``', 'x```')
                title([ char(filterStr(filterType)),' ',noisy,' Derivative Estimates, beta'])
                drawnow

                disp('PRESS ANY KEY')
                pause
            end
            kinematicDataX{noiseTest, filterType} = dx;
            kinematicDataY{noiseTest, filterType} = dy;


        end
    end

    kinematicResults = {};
    betaResults = [];
    VGFResults = [];

    for noiseTest = 2
        for filterType = 1:4;

            curvatureKinematic = [];
            curvatureGeometric = [];
            for regressTypes = 1:5
                % 1 regress with forced 0 intercept
                % 2 fitlm with forced 0 intercept
                % 3 fitlm with non zero intercept
                % 4 fitnlm Levenberg-Marquardt
                % 5 fitnlm Iteratively Reweighted Least Squares

                % clip the ends to avoid spurious edge elements from
                % differentiation and wide  Savitzky-Golay filters

                velocityX =  kinematicDataX{noiseTest, filterType}(edgeClip:end-edgeClip,2);
                velocityY =  kinematicDataY{noiseTest, filterType}(edgeClip:end-edgeClip,2);

                velocity = ( ( velocityX.^2 + velocityY.^2 ) .^0.5 );

                accelerationX =  kinematicDataX{noiseTest, filterType}(edgeClip:end-edgeClip,3);
                accelerationY =  kinematicDataY{noiseTest, filterType}(edgeClip:end-edgeClip,3);

                if displayCurvature % then we must calculate both for comparison

                    curvatureKinematic = curvatureKinematicEBR(velocityX, velocityY,accelerationX, accelerationY);

                    for ks = 1:length(accelerationY)
                        tripletXY = [ x(edgeClip+ks-1:edgeClip+ks+1) y(edgeClip+ks-1:edgeClip+ks+1) ]';
                        % we must take the absolute Menger Curvature to
                        % match the Kinematic Curvature
                        curvatureGeometric(ks,1) = abs(curvatureMengerEBR(tripletXY));
                    end

                    figure(99999)
                    clf
                    tiledlayout(3,1,'TileSpacing','Compact','Padding','Compact');
                    nexttile
                    plot(curvatureGeometric, 'r', 'DisplayName', 'Geometric Menger Curvature')
                    drawnow
                    legend
                    nexttile
                    plot(curvatureKinematic, 'b', 'DisplayName', 'Kinematic Curvature')
                    drawnow
                    legend
                    nexttile
                    plot(curvatureGeometric-curvatureKinematic,'k','DisplayName', 'Curvature Delta' )
                    legend
                    drawnow


                    if curvatureChoice == 1;
                        curvature = curvatureKinematic;
                    elseif curvatureChoice == 2;
                        curvature = curvatureGeometric;
                    end

                else
                    if curvatureChoice == 1
                        curvature = curvatureKinematicEBR(velocityX, velocityY,accelerationX, accelerationY);
                    else
                        for ks = 1:length(accelerationY)
                            tripletXY = [ x(edgeClip+ks-1:edgeClip+ks+1) y(edgeClip+ks-1:edgeClip+ks+1) ]';
                            % we must take the absolute Menger Curvature to
                            % match the Kinematic Curvature
                            curvature(ks,1) = abs(curvatureMengerEBR(tripletXY));
                        end
                    end
                end



                % cast them explicity for their role in the ensuing
                % regression
                responses = velocity;
                predictors = curvature;

                % we couch this in a try catch as sometimes iterative
                % regressions do not reach a result within the limited
                % iterations permitted.
                try
                    [beta, yGain, stats] = regressDataEBR(responses, predictors, regressTypes, LMSeeds, displayGraphs);
                    %   disp(['beta ' num2str(beta),' ', regressStr{regressTypes}, ' RESOLVES at NoiseSTD = ', num2str(noise(noises)), ' using ' filterStr{filterType}]);
                    betaResults(noiseTest, filterType, regressTypes) = beta;
                    VGFResults(noiseTest, filterType, regressTypes) = yGain;
                catch
                    disp([ regressStr{regressTypes}, ' DOESN`T resolve at NoiseSTD = ', num2str(noise(noises)), ' using ' filterStr{filterType}]);
                    betaResults(noiseTest, filterType, regressTypes) = NaN;
                    VGFResults(noiseTest, filterType, regressTypes) = NaN;
                end

                if displayGraphs
                    drawnow

                    disp('Regression Pause - Press Any Key');
                    %pause;
                end

            end
            %                             disp('PRESS ANY KEY')
            %                 pause
        end
    end

    resultDiff(noises,:) = betaResults(2, 1, :);
    resultBWDiff(noises,:) = betaResults(2, 2, :); % }
    resultDiffBW(noises,:) = betaResults(2, 3, :); % } the delta between these two is miniscule.
    resultSG(noises,:) = betaResults(2, 4, :);

    VGFDiff(noises,:) = VGFResults(2, 1, :);
    VGFBWDiff(noises,:) = VGFResults(2, 2, :);
    VGFDiffBW(noises,:) = VGFResults(2, 3, :);
    VGFSG(noises,:) = VGFResults(2, 4, :);

end


fig1 = figure(1);
set(groot,'defaultLineLineWidth',2)
fig1.WindowState = 'maximized';
clf
tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
sp(1) = nexttile;

colorsDot = ["r." "g." "b." "c." "m."];
colorsStr = ["r" "g" "b" "c" "m"];

for its = 1:5
    plot(noise,resultDiff(:,its), colorsDot(its), 'DisplayName',['Diff + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(resultDiff(:,its),length(resultDiff)/reps , reps)');
    plot(noise(1:length(resultDiff)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; %legend
ylabel('Beta'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise ',kinStr{curvatureChoice},' Diff'])
hold off

sp(2) = nexttile;
for its = 1:5

    plot(noise,resultBWDiff(:,its),colorsDot(its), 'DisplayName',['BWDiff + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(resultBWDiff(:,its),length(resultDiff)/reps , reps)');
    plot(noise(1:length(resultBWDiff)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; %legend
ylabel('Beta'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise ', kinStr{curvatureChoice},' BWDiff ', num2str(filterParamsBW)])
hold off

sp(3) = nexttile;
for its = 1:5

    plot(noise,resultDiffBW(:,its),colorsDot(its), 'DisplayName',['DiffBW + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(resultDiffBW(:,its),length(resultDiff)/reps , reps)');
    plot(noise(1:length(resultDiffBW)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; %legend
ylabel('Beta'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise','-',kinStr{curvatureChoice},' DiffBW ', num2str(filterParamsBW)])
hold off

sp(4) = nexttile;
for its = 1:5

    plot(noise,resultSG(:,its),colorsDot(its), 'DisplayName',['SG + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(resultSG(:,its),length(resultDiff)/reps , reps)');
    plot(noise(1:length(resultSG)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; legend
ylabel('Beta'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise','-',kinStr{curvatureChoice},' SG ', num2str(filterParamsSG)])
hold off

linkaxes(sp,'y');
ax = axis;
axis([ax(1:2) -1.1/3 0]);
figFilename1 = ['Beta_', num2str(max(maxVal)),'mm_',equiStr{equidistant},'-',kinStr{curvatureChoice},'-',noiseStr{noiseType} '-SG', num2str(filterParamsSG)];%'test.png'

fig2 = figure(2);
set(groot,'defaultLineLineWidth',2)
fig2.WindowState = 'maximized';
clf
tiledlayout(1,4,'TileSpacing','Compact','Padding','Compact');
sp(1) = nexttile;

colorsDot = ["r." "g." "b." "c." "m."];
colorsStr = ["r" "g" "b" "c" "m"];

for its = 1:5
    plot(noise,VGFDiff(:,its), colorsDot(its), 'DisplayName',['Diff + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(VGFDiff(:,its),length(VGFDiff)/reps , reps)');
    plot(noise(1:length(VGFDiff)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; %legend
ylabel('VGF'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise, Diff'])
hold off

sp(2) = nexttile;
for its = 1:5

    plot(noise,VGFBWDiff(:,its),colorsDot(its), 'DisplayName',['BWDiff + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(VGFBWDiff(:,its),length(VGFDiff)/reps , reps)');
    plot(noise(1:length(VGFBWDiff)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; %legend
ylabel('VGF'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise',' ',kinStr{curvatureChoice}, ' BWDiff ', num2str(filterParamsBW)])
hold off

sp(3) = nexttile;
for its = 1:5

    plot(noise,VGFDiffBW(:,its),colorsDot(its), 'DisplayName',['DiffBW + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(VGFDiffBW(:,its),length(VGFDiff)/reps , reps)');
    plot(noise(1:length(VGFDiffBW)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; %legend
ylabel('VGF'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise',' ',kinStr{curvatureChoice},' DiffBW ', num2str(filterParamsBW)])
hold off

sp(4) = nexttile;
for its = 1:5

    plot(noise,VGFSG(:,its),colorsDot(its), 'DisplayName',['SG + ', regressStr{its}])
    hold on
    resultAvg = nanmean(reshape(VGFSG(:,its),length(VGFDiff)/reps , reps)');
    plot(noise(1:length(VGFSG)/reps), resultAvg,colorsStr(its), 'DisplayName','')
end
grid on; legend
ylabel('VGF'); xlabel('StD Noise (mm)')
title([equiStr{equidistant},' + ',noiseStr{noiseType} ' noise,',' ',kinStr{curvatureChoice},' SG ', num2str(filterParamsSG)])
hold off

linkaxes(sp,'y');
ax = axis;
axis([ax(1:2) -1 7]);
figFilename2 = ['VGF_', num2str(max(maxVal)),'mm_',equiStr{equidistant},'-',kinStr{curvatureChoice},'-',noiseStr{noiseType} '-SG', num2str(filterParamsSG)];%'test.png'


if saveImage
folderSave = pwd;
folderSave = [folderSave(1:end-3) 'figures'];
    saveas(fig2,fullfile(folderSave,figFilename2),'png');
    saveas(fig1,fullfile(folderSave,figFilename1),'png');
end
