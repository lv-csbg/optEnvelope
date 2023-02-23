function [knockouts] = midsequentialOEReinserts(model, values, biomass, controlFlux1, desiredProduct, matchRev, K, ratio,MinProduct, timeLimit)

a=0;
MILP = struct([]);
var = [];
modelOld = model;
allActive = zeros(numel(model.rxns),length(values));
allModels = {};
for i = 1:length(values)
    a=a+1;
    model = modelOld;
    model=changeRxnBounds(model,biomass,controlFlux1(values(a)),'b');
    model=changeObjective(model,desiredProduct);
    s=optimizeCbModel(model);
    model=changeRxnBounds(model,desiredProduct,s.f*ratio,'b');
    tic
    [ActiveRxns, var, MILP]=minActiveRxns(model,matchRev,K, var, MILP, timeLimit);
    T = toc;
    test(1,i) = values(i);
    test(2,i) = T;
    test(3,i) = numel(ActiveRxns);
    if isempty(ActiveRxns)
        continue;
    end
    ids = findRxnIDs(model,ActiveRxns);
    allActive(1:length(ids),i) = ids;
    allActive(~any(allActive,2),:) = [];
    allActive(:,~any(allActive,1)) = [];
    allModels = [allModels;model];
end

allKnockouts = {};
[~,col] = size(allActive);
bestPoint = 0;
n=0;
for j = 1:col
    model = allModels{j};
    rxns = allActive(:,j);
    rxns(rxns==0)=[];
    rxns = ismember(model.rxns,[model.rxns(K);model.rxns(rxns)]);
    model.ub(rxns==0)=0;
    model=changeRxnBounds(model,desiredProduct,1000,'u');
    model=changeRxnBounds(model,desiredProduct,0,'l');
    model=changeRxnBounds(model,biomass,0,'b');
    s= optimizeCbModel(model,'min');MinB1=s.f;
    if MinB1>MinProduct
        model=changeRxnBounds(model,biomass,1000,'u');
        model=changeRxnBounds(model,biomass,0,'l');
        model2=model;n1=0;
        knockouts = {};
        
        % 3. Reduce the number of knockouts to minimum possible
        for i=1:numel(model.rxns)
            if ismember(i,find(rxns))==0
                model.lb(i)=modelOld.lb(i);
                model.ub(i)=modelOld.ub(i);
                solMin2 = optimizeCbModel(model,'min');
                if MinB1 - solMin2.f>10^-6
                    model.ub(i)=model2.ub(i);model.lb(i)=model2.lb(i);
                    n1=n1+1;knockouts(n1,1)=model.rxns(i);
                end
            end
        end
        n=n+1;
        allKnockouts{n} = knockouts;
        if isempty(knockouts) 
            continue;
        elseif bestPoint == 0
            bestPoint = n;
        elseif length(allKnockouts{bestPoint})>length(knockouts)
            bestPoint = n;
        end
    end
end
if bestPoint==0
    knockouts = {};
else
    knockouts = allKnockouts{bestPoint};
end