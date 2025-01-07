% Demo Code for One-Third Power Law Principled Protocol
%
% Velocity Gain Factor and Beta Power Law exponents are then calculated via
% choices of differentiation, filtering, and regression
% as outlined in Exp Brain Res review paper of Fraser et al., 2024
%
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
% requires functions subfolder with bespoke;
% - curvatureKinematicEBR.m
% - differentiateKinematicsEBR.m
% - regressDataEBR.m
% - importfileZarandi.m
%
% and data files from Zarandi et al. 2023 DOI: 10.1038/s41598-023-34861-x
% https://github.com/lucaoneto/IJCNN2022_Ellipses/tree/main/data
%
% v0.1.0
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 


clear all
close all

addpath(genpath('functions'))
addpath(genpath('req'))

subjectsFiles =     { 'subject0.txt' 'subject1.txt'   'subject2.txt'  'subject3.txt'   'subject4.txt'  'subject5.txt'  ...
    'subject6.txt'  'subject7.txt'   'subject8.txt'   'subject9.txt'  'subject10.txt' ...
    'subject11.txt'  'subject12.txt' 'subject13.txt' };

rng("default") % make it repeatable
regenData = 1; % 1 re-create the results from raw data
saveImage = 1; % create publication image
displayGraphs = 0; % display intermediate graphs
coherenceCheck = 0; % display a nice time based graph of the data being drawn to check for anti clockwise

if regenData

    coherenceCheck = 0;
    pauseCheck = 0;
    dataBig = [];
    localDir = pwd;
    masterDir = localDir(1:end-3);
    edgeClip = 20;

    LMseeds = [1 -1/3];
    displayGraphs = 0 ;
    limitBreak = 0;

    beta = NaN(14,10,4,5);
    VGF = NaN(14,10,4,5);

    for its = 1:14

        [subjID, dominant_hand, speed, trial, nellipse, x, y, Velocity, Pressure, Radius, t]...
            = importfileZarandi(fullfile( [masterDir '/data/_RawDataZarandi/IJCNN2022_Ellipses-main/IJCNN2022_Ellipses-main/data/' subjectsFiles{its}]));

        subjID(1);

        % identify trials with only the dominant hand... as the ND will need to be
        % reversed for anti clockwise and we just want fast Beta results 

        dominantHand = find(dominant_hand);
        if displayGraphs 
                clf
                plot(trial(dominantHand)); title('Slow - Natural - Fast')
                drawnow
        end
        %between the second and thirds groups of 00s

        % now find the trial 0s as this segments the sets of trials
        % the middle ones the natural frequency trials

        % https://uk.mathworks.com/matlabcentral/answers/181400-how-can-i-find-consecutive-ones-and-the-number-of-ones-in-each-block
        % x=[5 1 1 1 1 1 7 1 1 5];
        % start1 = strfind([0,x==1],[0 1]);
        % end1 = strfind([x==1,0],[1 0]);
        % end1 - start1 + 1
        % ans =
        %      5     2

        %    to avolid entering the 3rd format
        start1Max3 = strfind([0,trial(dominantHand)'==0],[0 1]);
        start1Max = start1Max3(3);

        for trialNum = 0:9
            disp(['Participant ', num2str(its), '- Trial ', num2str(trialNum)])
            start1 = strfind([0,trial(dominantHand)'==trialNum],[0 1]);
            end1 = strfind([trial(dominantHand)'==trialNum,0],[1 0]);

            xDominant = x(dominantHand);
            yDominant = y(dominantHand);
            tDominant = t(dominantHand);

if displayGraphs
            figure(10)
            subplot(3,1,1)
            plot(trial(dominantHand)); title('Slow - Natural - Fast(?)')
            hold on
end

            if (length(start1)>=2)  && (start1(2) < start1Max)

                xTrial = xDominant(start1(2):end1(2));
                yTrial = yDominant(start1(2):end1(2));
                tTrial = tDominant(start1(2):end1(2));
if displayGraphs 
                figure(10)
                subplot(3,1,2)
                plot(xTrial, yTrial)

                subplot(3,1,1)
                plot(trial(dominantHand)); title('Slow - Natural - Fast(?)')
                hold on
                scatter(start1(2):100:end1(2), trialNum);
                hold off
                drawnow;
end
                if coherenceCheck %% draw the data as if were beign gathered live
                    figure(99)
                    for tCount = 1:10:length(tTrial)
                        scatter(xTrial(tCount), yTrial(tCount))
                        axis([0   30.0000  -30.0000  0.0000])
                        drawnow
                    end
                    disp('Press Any Key')
                    if pauseCheck
                        pause
                    end

                end

                %% we can parcel this data off when we save

                dataBig{its, trialNum+1}.xTraj = xTrial;
                dataBig{its, trialNum+1}.yTraj = yTrial;
                dataBig{its, trialNum+1}.tTime = tTrial;
                fs = 100;



                for diffType = 1:3

                    if diffType == 1 % Zarandi Filter
                        filterType = 2; %Low pass BW
                        filterParams = [1 0.07 1];% order Fc zerolag switch ON
                    elseif diffType == 2 % Zarandi Filter Corrected?
                        filterType = 2; %Low pass BW
                        filterParams = [2 10 1];% order Fc zerolag switch ON
                    elseif diffType == 3 % Zarandi Filter
                        filterType = 4; %Savitzky-Golay
                        filterParams = [4 17];% order framelength a la Crenna et al,
                    end

                    [dx, dy] = differentiateKinematicsEBR(xTrial, yTrial, filterType, filterParams, fs);

                    velocityX =  dx(edgeClip:end-edgeClip,2);
                    velocityY =  dy(edgeClip:end-edgeClip,2);
                    velocity = ( ( velocityX.^2 + velocityY.^2 ) .^0.5 );

                    accelerationX =  dx(edgeClip:end-edgeClip,3);
                    accelerationY =  dy(edgeClip:end-edgeClip,3);

                    curvatureK = curvatureKinematicEBR(velocityX, velocityY,accelerationX, accelerationY);

                    responses = velocity;
                    predictors = curvatureK;
if displayGraphs 
                    figure(10)
                    subplot(3,1,3)
                    plot(curvatureK)
                    hold on
                    plot(velocity)
                    drawnow
                    hold off
end

                    for regressType = 3:5
                        % 3 fitlm with non zero intercept
                        % 4 fitnlm Levenberg-Marquardt
                        % 5 fitnlm Iteratively Reweighted Least Squares

                        betaLocal = NaN;
                        yGain = NaN;

                        [betaLocal,yGain, stats] = regressDataEBR(responses, predictors, regressType, LMseeds, displayGraphs, limitBreak);
if displayGraphs 
                        title([' Beta = ',num2str(betaLocal),' Filter= ',num2str(diffType),' Regress = ',num2str(regressType)])
end
                        disp([' Beta = ',num2str(betaLocal),' Filter= ',num2str(diffType),' Regress = ',num2str(regressType)])

                        %% use the Linear Regression to update the seed for the non-linear regression
                        if regressType == 3
                            LMseeds = [yGain betaLocal];
                        end

                        beta(its, trialNum+1, diffType, regressType) = betaLocal;
                        VGF(its, trialNum+1, diffType, regressType) = yGain;
                    end
                end
            end
        end
    end

    save % default matlab.mat

else
    load % default matlab.mat
end
% we now have beta and VGF
%  1-14 pIDs  0 upto 9 trial, 1-3 filtering and diff, 3-5 fitlm, LMLS, IRLS

%we want to see if the Beta differs between different methods

for diffType  =2:3
    for regressType = 3:5

        betaArray = beta(:,:,diffType,regressType);
        % take the subset of the data for all pIDs and trials, for one set of filter and regress options
        betaSerial =  betaArray(:); % flatten it to take the mean and std nicely

        betaSerialCell{diffType, regressType} = betaSerial;

        betaMean(diffType, regressType) = nanmean(betaSerial);
        betaStd(diffType, regressType) = nanstd(betaSerial);

    end
end

% https://stackoverflow.com/questions/47502438/matlab-combination-boxplot-of-different-length-vectors

g1 = ones(size(betaSerialCell{2,3})) * 1;
g2 = ones(size(betaSerialCell{2,4})) * 2;
g3 = ones(size(betaSerialCell{2,5})) * 3;
g4 = ones(size(betaSerialCell{3,3})) * 4;
g5 = ones(size(betaSerialCell{3,4})) * 5;
g6 = ones(size(betaSerialCell{3,5})) * 6;

figure(11)
set(groot,'defaultLineLineWidth',1.0)
boxplot([betaSerialCell{2,3}; betaSerialCell{2,4}; betaSerialCell{2,5}; ...
    betaSerialCell{3,3}; betaSerialCell{3,4}; betaSerialCell{3,5}],  ...
    [g1;g2;g3;g4;g5;g6], 'Notch','on', 'Labels', ...
    {'BW LR', 'BW LMLS', 'BW IRLS', 'SG LR', 'SG LMLS', 'SG IRLS'})
grid on
yl = yline([-1/3],'-.r', '-1/3', 'LineWidth',3)
yl.LabelHorizontalAlignment = 'left';
title({'Comparison of variant Filtering, Differentiation and Regression.', 'Data from Zarandi et al. 2023'})
ylabel('Empirical Beta')
xlabel({ ' ','Filter Types; BW Butterworth, SG Savitzky-Golay',...
    'Regression Types; LR Linear, LMLS Levenberg-Marquardt, IRLS Iteratively Reweighted ', ' '})

% https://uk.mathworks.com/matlabcentral/answers/175193-creating-sigstar-in-bar-graph

 fig11 = gcf;

    figFilename11 = ['Zarandi','_paperFigure'];
    if saveImage
        folderSave = pwd;
        folderSave = [folderSave(1:end-3) 'figures'];
        % saveas(fig11,fullfile(folderSave,figFilename11),'png');
        % Requires R2020a or later
        plotedit(fig11,'on');
        disp('PRESS ANY KEY ONCE PLOT ELEMENTS MOVED TO FINAL POSITION!')
        pause
        plotedit(fig11,'off');
        exportgraphics(fig11,fullfile(folderSave,[figFilename11,'.png']),'Resolution',600);
    end
    close all
