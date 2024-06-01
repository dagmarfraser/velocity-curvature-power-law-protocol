function [beta,yGain, stats] = regressDataEBR(responses, predictors, regressType, LMseeds, displayGraphs, limitBreak)
%% demo code for Exp Brain Res review paper of Fraser et al., 2024
% Created May 2024
% Correspondence Dagmar Scott Fraser
% d.s.fraser@bham.ac.uk
%
%% inputs
% regress reponses / predictors
% standard MATLAB is regress(y,X) where y is 'responses'
% and X is the 'predictor'
% We implement several 'regresstype' with common function inputs
% 1 regress with forced 0 intercept
% 2 fitlm with forced 0 intercept
% 3 fitlm with non zero intercept
% 4 fitnlm Levenberg-Marquardt
% 5 fitnlm Iteratively Reweighted Least Squares
%
% for model fitting such as 4 & 5 we require LMseeds
% e.g. [VGF Beta] = [1 -1/3]
%
% displayGraphs - 0 - no report, 1 - show regression
%
% limitBreak - if limitBreak = 1 then we keep the last value of the
% Nonlinear fit... otherwise we return Nan if it has not settled within
% MaxIter = 100
%
%% outputs
% estimates of Beta and VGF (not log projected) , associated stats

logy = log(responses); % log data - log(v)
logx = log(predictors); % log data - log(k)

y = responses; % (v)
x = predictors; % (k)

stats = [];

switch regressType
    case 1 % zero intercept regression

        [b,bint,r] = regress(logy,logx); % regress(y,X)
        beta = b;
        yGain = NaN;
        stats = { bint, r };

    case 2 % should give identical results to case 0
        %mdl = fitlm(X,y) - note they swap the order of X and y

        % 'standard' regression y = ax
        mdlFbeta = fitlm(logx,logy, 'Intercept',false);
        beta = mdlFbeta.Coefficients.Estimate;
        yGain = NaN;

        if displayGraphs
            figure(1001)
            clf
            plot(mdlFbeta)
            title(['Regression FITLM (-intercept) log v = ',num2str(mdlFbeta.Coefficients.Estimate), ' * log k'])
            xlabel('log k')
            ylabel('log v')
        end

    case 3
        % stardard regression with potentialo for non zero intercept i.e. y = a + bx
        mdlTbeta = fitlm(logx,logy, 'Intercept',true);
        %beta = [mdlTbeta.Coefficients.Estimate(1) mdlTbeta.Coefficients.Estimate(2)];
        beta = mdlTbeta.Coefficients.Estimate(2);
        yGain = exp(mdlTbeta.Coefficients.Estimate(1));
        if displayGraphs
            figure(1002)
            clf
            plot(mdlTbeta)
            title(['Regression FITLM (+intercept) log v = ',num2str(mdlTbeta.Coefficients.Estimate(1)), ' + ',num2str(mdlTbeta.Coefficients.Estimate(2)),' * log k'])
            xlabel('log k')
            ylabel('log v')
            hold off
        end


    case 4
        % For nonrobust estimation, nlinfit uses the Levenberg-Marquardt
        % nonlinear least squares algorithm.. nonrobust is the default, and
        % only 100 iterations.
        % https://uk.mathworks.com/help/stats/nlinfit.html

        beta0 = LMseeds(2); % e.g. -1/3 for the ellipse
        yGain0 = LMseeds(1); % anecdotally <1, so after log it will be <0

        modelFunc = @(beta, s)( beta(1)*s.^(beta(2)));
        % speed = yGain * curvature ^-Beta
        % NB Beta return as negative
        % to match all the other mdl / regressions

        %% we want to not return guesses - Iteration limit exceeded.  Returning results from final iteration.
        lastwarn('')

        initialGuess = [yGain0 ; beta0]; % initial guesses
        %[betaNLF,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y, modelFunc, initialGuess)
        % use the below version which returns the model.
        mdlLMbeta = fitnlm(x,y, modelFunc, initialGuess);
        % this may return Inf and NaN which will break

        if isempty(lastwarn) || limitBreak % we have resolved within MaxIter OR we don't care
            if limitBreak
                disp('Returning Best Guess for Beta and VGF')            %beta = [log(mdlLMbeta.Coefficients.Estimate(1)) mdlLMbeta.Coefficients.Estimate(2)];
            end
            %beta = [log(mdlLMbeta.Coefficients.Estimate(1)) mdlLMbeta.Coefficients.Estimate(2)];
            beta = mdlLMbeta.Coefficients.Estimate(2);
            yGain = (mdlLMbeta.Coefficients.Estimate(1));
        else
            % 'Iteration limit exceeded.  Returning results from final iteration.'
            disp('Returning NaN for Beta and VGF')
            beta = NaN;
            yGain = NaN;
        end
        % to match the log log regression above we need to log the yGain.

        if displayGraphs
            figure(1003)
            clf
            pars = mdlLMbeta.Coefficients.Estimate;
            yMdl = modelFunc(pars,x);
            subplot(2,1,1)

            scatter(x,y, '.')
            hold on
            plot(x,yMdl, '.');
            title(['Target beta0= ',num2str(beta0),' LevenbergMarquardt FITNLM with model func y = a * x^b i.e. in non log log'])
            xlabel('k')
            ylabel('v')
            subplot(2,1,2)
            scatter(logx,logy, 'blue', 'x')
            hold on
            scatter(logx,log(yMdl), 'red', '+', 'LineWidth',0.25)
            %            title(['LevenbergMarquardt FITNLM (+intercept) v = ',num2str(mdlLMbeta.Coefficients.Estimate(1)), ' * k ^ ',num2str(mdlLMbeta.Coefficients.Estimate(2))])
            title(['LevenbergMarquardt FITNLM (+intercept) log v = ',num2str(log(mdlLMbeta.Coefficients.Estimate(1))), ' + ',num2str(mdlLMbeta.Coefficients.Estimate(2)), ' * log k'])
            xlabel('log k')
            ylabel('log v')
            hold off
        end

    case 5
        % For robust estimation, nlinfit uses the Iteratively Reweighted Least Squares
        % nonlinear least squares algorithm.. nonrobust is the default, and
        % only 100 iterations.

        beta0 = LMseeds(2); % e.g. -1/3 for the ellipse
        yGain0 = LMseeds(1); % anecdotally <1, so after log it will be <0

        modelFunc = @(beta, s)( beta(1)*s.^(beta(2)));
        % speed = yGain * curvature ^-Beta
        % note we let Beta return as negative
        % toi match all the other mdl / regressions

        %% we want to not return guesses - Iteration limit exceeded.  Returning results from final iteration.
        lastwarn('')

        initialGuess = [yGain0 ; beta0]; % initial guesses
        %[betaNLF,R,J,CovB,MSE,ErrorModelInfo] = nlinfit(x,y, modelFunc, initialGuess)
        % use the below version which returns the model.
        % Set robust fitting options.

        opts = statset('nlinfit');
        opts.RobustWgtFun = 'bisquare';

        mdlLMbeta = fitnlm(x,y, modelFunc, initialGuess, 'Options', opts);

        if isempty(lastwarn) || limitBreak % we have resolved within MaxIter OR we don't care
            if limitBreak
                disp('Returning Best Guess for Beta and VGF')            %beta = [log(mdlLMbeta.Coefficients.Estimate(1)) mdlLMbeta.Coefficients.Estimate(2)];
            end
            beta = mdlLMbeta.Coefficients.Estimate(2);
            yGain = (mdlLMbeta.Coefficients.Estimate(1));
            % to match the log log regression above we need to log the yGain.
        else
            % lastwarn == 'Iteration limit exceeded.  Returning results from final iteration.'
            disp('Returning NaN for Beta and VGF')
            beta = NaN;
            yGain = NaN;
        end

        if displayGraphs
            figure(1004)
            clf
            pars = mdlLMbeta.Coefficients.Estimate;
            yMdl = modelFunc(pars,x);
            subplot(2,1,1)

            scatter(x,y, '.')
            hold on
            plot(x,yMdl, '.');
            title(['Target beta0= ',num2str(beta0),' Iteratively Reweighted Least Squares FITNLM with model func y = a * x^b i.e. in non log log'])
            xlabel('k')
            ylabel('v')
            subplot(2,1,2)
            scatter(logx,logy, 'blue', 'x')
            hold on
            scatter(logx,log(yMdl), 'red', '+', 'LineWidth',0.25)
            title(['Iteratively Reweighted Least Squares FITNLM (+intercept) log v = ',num2str(mdlLMbeta.Coefficients.Estimate(1)), ' + ',num2str(mdlLMbeta.Coefficients.Estimate(2)),' * log k'])
            xlabel('log k')
            ylabel('log v')
            hold off
        end


end