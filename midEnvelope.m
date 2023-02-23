function [midknockouts,ratioF] = midEnvelope(model, values, biomass, controlFlux1, desiredProduct, matchRev, K,knockouts,nMidpoints,~, timeLimit)
model1=model;
rxns=ismember(model.rxns,knockouts);
maxKO = numel(knockouts) + 2;
model.ub(rxns)=0;
model.lb(rxns)=0;
model=changeObjective(model,biomass);
s= optimizeCbModel(model,'max');MaxBiomass=s.f;
model=changeRxnBounds(model,biomass,MaxBiomass,'b');
model=changeObjective(model,desiredProduct);
s= optimizeCbModel(model,'max');midProduct=s.f;
model=changeRxnBounds(model,biomass,0,'b');
s= optimizeCbModel(model,'min');MinProduct=s.f;
s= optimizeCbModel(model,'max');MaxProduct=s.f;
minRatio=midProduct/MaxProduct;
ratio = linspace(minRatio,1,nMidpoints);n=0;
for i=1:numel(ratio)-1
    %if isempty(numKO)
    [knockouts1] = midsequentialOEReinserts(model1, values, biomass, controlFlux1, desiredProduct, matchRev, K, ratio(i),MinProduct, timeLimit);
    %else
    %    [knockouts1] = midOEReinserts(model1, values, biomass, controlFlux1, desiredProduct, matchRev, K, ratio(i), numKO, timeLimit);
    %end
    if ~isempty(knockouts1)
        if numel(knockouts1)<=maxKO
            n=n+1;
            midknockouts{n}=knockouts1;
            ratioF(n)=ratio(i);
            number(n)=numel(knockouts1);
        end
    end
end
midknockouts{n+1}=knockouts;number(n+1)=numel(knockouts);ratioF(n+1)=1;
for i=1:n+1
    model=model1;
    rxns=ismember(model.rxns,midknockouts{i});
    model.ub(rxns)=0;
    model.lb(rxns)=0;
    model=changeObjective(model,biomass);
    s= optimizeCbModel(model,'max');MaxBiomass(i)=s.f;
    model=changeRxnBounds(model,biomass,MaxBiomass(i),'b');
    model=changeObjective(model,desiredProduct);
    s= optimizeCbModel(model,'max');midProduct(i)=s.f;
    model=changeRxnBounds(model,biomass,0,'b');
    s= optimizeCbModel(model,'min');MinProduct(i)=s.f;
    s= optimizeCbModel(model,'max');MaxProduct(i)=s.f;
end
T = table(ratioF',number',MaxBiomass',midProduct',MaxProduct',MinProduct','VariableNames',{'Ratio','No. of KO','Growth rate at MPP','Prod. rate at MPP','Max prod.','Min prod.'}) 


