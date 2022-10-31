function [knockouts, ActiveRxns] = optEnvelope(model,biomass,desiredProduct,desiredProductName,matchRev,K,prodMol)
%optEnvelope finds minimum active reactions and uses a sequential method to
% reinsert reactions and propose minimum knockouts supporting optimal
% envelope
%
%INPUT
% model              COBRA model structure
% biomass            Biomass reaction name
% desiredProduct     Reaction name of desired product
% desiredProductName Desired product name
% matchRev           Matching of forward and backward reactions of a reversible
%                    reaction
% K                  List of reactions that cannot be selected for knockout
% prodMol            Molar mass of product
%
%OUTPUT
% knockouts          List of which reactions to support the optimal
% ActiveRxns         List of all active reactions
%
%Notes
% It should be mentioned that a figure (desired product versus biomass)
% including plots for wild-type and opt enveople is presented after
% running optEnvelope
%
% created by  Ehsan Motamedian       09/02/2022
% modified by Kristaps Berzins       31/10/2022
model1=model;

% 1. create wild-type envelope
solMin = optimizeCbModel(model,'min');
solMax = optimizeCbModel(model,'max');
controlFlux1 = linspace(solMin.f,solMax.f,100)';
model=changeObjective(model,desiredProduct);
for i=1:numel(controlFlux1)
    model=changeRxnBounds(model,biomass,controlFlux1(i),'b');
    s= optimizeCbModel(model,'min');Min1(i,1)=s.f;
    s= optimizeCbModel(model,'max');Max1(i,1)=s.f;
end
figure()
hold on
p1=plot(controlFlux1/10*1000/180,Max1/180*prodMol/10,'b','LineWidth',2);
plot(controlFlux1/10*1000/180,Min1/180*prodMol/10,'b','LineWidth',2)
xlabel('Biomass yield (g/g)')
ylabel([desiredProductName,' production yield (g/g)'])

% 2. Calculate minimum active reactions
model=changeRxnBounds(model,biomass,controlFlux1(8),'b');
s=optimizeCbModel(model);
model=changeRxnBounds(model,desiredProduct,s.f,'b');
ActiveRxns=minActiveRxns(model,matchRev,K);
rxns=ismember(model.rxns,[model.rxns(K);ActiveRxns]);
model.ub(rxns==0)=0;
model=changeRxnBounds(model,desiredProduct,1000,'u');
model=changeRxnBounds(model,desiredProduct,0,'l');
model=changeRxnBounds(model,biomass,0,'b');
s= optimizeCbModel(model,'min');MinB1=s.f;
model=changeRxnBounds(model,biomass,1000,'u');
model=changeRxnBounds(model,biomass,0,'l');
model2=model;n1=0;

% 3. Reduce the number of knockouts to minimum possible
for i=1:numel(model.rxns)
    if ismember(i,find(rxns))==0
        model.lb(i)=model1.lb(i);
        model.ub(i)=model1.ub(i);
        solMin2 = optimizeCbModel(model,'min');
        if MinB1 - solMin2.f>10^-6
            model.ub(i)=model2.ub(i);model.lb(i)=model2.lb(i);
            n1=n1+1;knockouts(n1,1)=model.rxns(i);
        end
    end
end

% 4. Plot minimal production envelope
model=changeObjective(model,biomass);
prID =findRxnIDs(model,desiredProduct);
solMin = optimizeCbModel(model,'min');
solMax = optimizeCbModel(model,'max');
controlFlux = linspace(solMin.f,solMax.f,100)';
model=changeObjective(model,desiredProduct);

for i=1:numel(controlFlux)
    model=changeRxnBounds(model,biomass,controlFlux(i),'b');
    s= optimizeCbModel(model,'min');Min(i,1)=s.f;
    s= optimizeCbModel(model,'max');Max(i,1)=s.f;
end
Min(end) = solMax.x(prID);
Max(end) = solMax.x(prID);
p3=plot(controlFlux/180*1000/10,Max/180*prodMol/10,'r','LineWidth',2);
plot(controlFlux/180*1000/10,Min/180*prodMol/10,'r','LineWidth',2)
legend([p1 p3],{'Wild-type','optEnvelope'})
hold off
