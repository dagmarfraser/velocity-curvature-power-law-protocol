% PowerLaw

load

% beta(pIDs, trialID, filtdiff 2 or 3, regress 3 4 or 5
% convert to table

tableT = [];
for pIDs = 1: 14
    for trials = 1:10
        for filtDiffs = 2:3
            for regresses = 3:5
                if isnan(beta(pIDs, trials, filtDiffs, regresses))
                tableT = [tableT ;NaN NaN NaN NaN NaN];               
                else
                tableT(end+1, 1) = beta(pIDs, trials, filtDiffs, regresses);
                tableT(end, 2) = pIDs;
                tableT(end, 3) = trials;
                tableT(end, 4) = filtDiffs;
                tableT(end, 5) = regresses;
                end
            end
        end
    end
end

tableTrue = array2table(tableT,...
    'VariableNames',{'Beta','ID','Trial', 'FiltDiff','Regress'});

tableTrue.ID = nominal(tableTrue.ID);
tableTrue.Trial = nominal(tableTrue.Trial);
tableTrue.FiltDiff = nominal(tableTrue.FiltDiff);
tableTrue.Regress = nominal(tableTrue.Regress);

lme = fitlme(tableTrue, 'Beta~FiltDiff*Regress+(1|ID) +(1|Trial')

