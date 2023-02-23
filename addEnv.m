function line = addEnv(origModel,biomass, desiredProduct, KnockOuts, colour, prodMol, subUptake, molarSum)
% addEnv adds envelope to figure
%
%INPUT
% origModel      COBRA model structure
% biomass        Biomass reaction name
% desiredProduct Reaction name of desired product
% KnockOuts      List of knockouts for production envelope
% colour         Short name for colour of line to plot (default: red)
% prodMol        Molar mass of target product for yield plot
% subUptake      Uptake of substrate
% molarSum       Molarmass of substrate
%
%OUTPUT
% line           Line data for plot function
%
%NOTES
% Sometimes last point of envelope drops to zero (might be rounding error)
% but this function connects last points of lines so the graph creates
% continuous line.
%
%Created by Kristaps Berzins    31/10/2022

if nargin < 4
    KnockOuts = {};
end

if nargin < 5
    colour = 'r';
end

if nargin < 6
    prodMolB = false;
    subUptake = [];
else
    prodMolB = true;
end

if nargin < 7
    subUptake = 10;
end

if nargin < 8
    molarSum = 180;
end

model = origModel;

rxns=ismember(model.rxns,KnockOuts);
model.ub(rxns)=0;
model.lb(rxns)=0;

solMin = optimizeCbModel(model,'min');
solMax = optimizeCbModel(model,'max');
controlFlux1 = linspace(solMin.f,solMax.f,100)';
if nnz(controlFlux1) == 0
    return;
end
model=changeObjective(model,desiredProduct);
%prodMax = optimizeCbModel(model,'max');
for i=1:numel(controlFlux1)
    model=changeRxnBounds(model,biomass,controlFlux1(i),'b');
    s= optimizeCbModel(model,'min');Min1(i,1)=s.f;
    s= optimizeCbModel(model,'max');Max1(i,1)=s.f;
end
if Max1(end) > Min1(end)
    Min1(end) = Max1(end);
elseif Max1(end) < Min1(end)
    Max1(end) = Min1(end);
elseif abs(Max1(end) - Min1(end)) < 1e-6
    Max1(end) = Min1(end-1);
    Min1(end) = Min1(end-1);
end

if prodMolB
    line=plot(controlFlux1/subUptake*1000/molarSum,Max1/molarSum*prodMol/subUptake,'color',colour,'LineWidth',2);
    plot(controlFlux1/subUptake*1000/molarSum,Min1/molarSum*prodMol/subUptake,'color',colour,'LineWidth',2)
else
    line=plot(controlFlux1,Max1,'color',colour,'LineWidth',2);
    plot(controlFlux1,Min1,'color',colour,'LineWidth',2)
end

