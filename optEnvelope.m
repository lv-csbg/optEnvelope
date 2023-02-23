function [knockouts, midknockouts] = optEnvelope(model, desiredProduct, varargin)
% optEnvelope finds minimum active reactions and uses a MILP method to
% reinsert reactions and propose minimum knockouts supporting optimal
% envelope
%
%   EXAMPLE: [knockouts, midKOs] = optEnvelope(model, 'EX_ac_e','timeLimit', 600, 'protectedRxns', {'H2Ot_f','H2Ot_b'}, 'nMidP', 15);
%
% INPUT
%  model              COBRA model structure
%  desiredProduct     Reaction ID of desired product
%
%  protectedRxns      (opt) Aditional reactions to ignore (must be in irreversible form) (default: {})
%  numKO              (opt) Number of reactions to remove for final result (default: [])
%  prodMol            (opt) Molar mass of product for yield plot (g/mol) (default: [])
%  values             (opt) Number of points to check along the edge for best envelope  (default: 10)
%  nMidP              (opt) Number of mid points between most probable
%                           optimal point and maximum product (default: 5)
%  midPointsOn        (opt) Baniry value if algorithm should calculate
%                           midpoints (default: true)
%  timeLimit          (opt) Time limit for gurobi optimization (default: inf)
%  printLevel         (opt) Print level for gurobi optimization (default: 0)
%  drawEnvelope       (opt) Binary value if algorithm should draw envelopes
%                           (default: true)
%
% OUTPUT
%  knockouts          List of which reactions to remove for the optimal
%                     envelope
%  midknockouts       List of reactions to remove for midpoint envelopes
%
% NOTES
%  It should be mentioned that a figure (desired product versus biomass)
%  including plots for wild-type and opt enveople is presented after
%  running optEnvelope
%
% created by  Ehsan Motamedian     02/09/2022
% modified by Kristaps Berzins     06/12/2022
% modified by Ehsan Motamedian     25/01/2023 switch to middle points was added


%% 0. Set parameters
parser = inputParser();
parser.addRequired('model', @(x) isstruct(x) && isfield(x, 'S') && isfield(model, 'rxns')...
    && isfield(model, 'mets') && isfield(model, 'lb') && isfield(model, 'ub') && isfield(model, 'b')...
    && isfield(model, 'c'))
parser.addRequired('desiredProduct', @(x) ischar(x))
parser.addParameter('protectedRxns', {}, @(x) iscell(x) && ismatrix(x));
parser.addParameter('numKO', [], @(x) isnumeric(x));
parser.addParameter('prodMol', [], @(x) isnumeric(x));
parser.addParameter('values', 10, @(x) isnumeric(x));
parser.addParameter('nMidP', 5, @(x) isnumeric(x));
parser.addParameter('midPointsOn', true, @(x) islogical(x));
parser.addParameter('timeLimit', inf, @(x) isnumeric(x));
parser.addParameter('printLevel', 0, @(x) isnumeric(x) || islogical(x));
parser.addParameter('drawEnvelope', true, @(x) islogical(x));

parser.parse(model, desiredProduct, varargin{:});
modelOri = parser.Results.model;
desiredProduct = parser.Results.desiredProduct;
protectedRxns = parser.Results.protectedRxns;
numKO = parser.Results.numKO;
prodMol = parser.Results.prodMol;
values = parser.Results.values;
nMidP = parser.Results.nMidP;
midPointsOn = parser.Results.midPointsOn;
timeLimit = parser.Results.timeLimit;
printLevel = parser.Results.printLevel;
drawEnvelope = parser.Results.drawEnvelope;

if isempty(prodMol)
    prodMolB = false;
else
    prodMolB = true;
end

values = (3:1:values);
[model,matchRev,~,~] = convertToIrreversible(modelOri);
K = findExcRxns(model);K = model.rxns(K);K = findRxnIDs(model,K);
if ~isempty(protectedRxns)
    KOid = findRxnIDs(model,protectedRxns);
    if any(KOid == 0)
        disp('At least one of reactions are not in the model - ignoring those')
        KOid(KOid==0) = [];
    end
    K = [K;KOid'];
    K = unique(K);
end
model1=model;

biomass = model.rxns(model.c == 1);
biomass = biomass{1};

desiredProductName = model.metNames(logical(abs(model.S(:, findRxnIDs(model, desiredProduct)))));
desiredProductName = desiredProductName{1};

if prodMolB
    input = model.rxns(model.ub < max(model.ub));
    numSub = size(input,1);
    if numSub > 1
        prompt = {'Choose substrate reaction'};
        answer = listdlg('PromptString',prompt,'SelectionMode','single','ListString',input);
    else
        answer = 1;
    end
    subUptake = model.ub(findRxnIDs(model, input(answer)));
    formula = model.metFormulas(logical(abs(model.S(:, findRxnIDs(model, input(answer))))));
    formula = formula{:};
    C = 12;
    H = 1;
    O = 16;
    indC = strfind(formula, 'C');
    indH = strfind(formula, 'H');
    indO = strfind(formula, 'O');
    C = C*str2double(formula(indC+1:indH-1));
    H = H*str2double(formula(indH+1:indO-1));
    O = O*str2double(formula(indO+1:end));
    molarSum = C + H + O;
end

%% 1. Create wild-type envelope
solMin = optimizeCbModel(model,'min');
solMax = optimizeCbModel(model,'max');
controlFlux1 = linspace(solMin.f,solMax.f,100)';
model=changeObjective(model,desiredProduct);
for i=1:numel(controlFlux1)
    model=changeRxnBounds(model,biomass,controlFlux1(i),'b');
    s= optimizeCbModel(model,'min');Min1(i,1)=s.f;
    s= optimizeCbModel(model,'max');Max1(i,1)=s.f;
end

if drawEnvelope
    figure()
    hold on

    if prodMolB
        p1=plot(controlFlux1/subUptake*1000/molarSum,Max1/molarSum*prodMol/subUptake,'b','LineWidth',2); %modify here for substrate
        plot(controlFlux1/subUptake*1000/molarSum,Min1/molarSum*prodMol/subUptake,'b','LineWidth',2)
        xlabel('Biomass yield (g/g)')
        ylabel([desiredProductName,' production yield (g/g)'])
    else
        p1=plot(controlFlux1,Max1,'b','LineWidth',2);
        plot(controlFlux1,Min1,'b','LineWidth',2)
        xlabel('Biomass(1/h)')
        ylabel([desiredProductName,' production (mmol/gDCW/h)'])
    end
end

%% 2. Reduce the number of knockouts to minimum possible and calculate midEnvelopes
warning off
if isempty(numKO)
    [knockouts] = sequentialOEReinserts(model, values, biomass, controlFlux1, desiredProduct, matchRev, K, model1, timeLimit);
    if midPointsOn
    	[midknockouts,ratio] = midEnvelope(model1, values, biomass, controlFlux1, desiredProduct, matchRev, K,knockouts,nMidP, numKO, timeLimit);
    end
else
    [knockouts] = findOEReinserts(model1, values, biomass, controlFlux1, desiredProduct, matchRev, K, numKO, timeLimit, printLevel);
    if midPointsOn
    	[midknockouts,ratio] = midEnvelope(model1, values, biomass, controlFlux1, desiredProduct, matchRev, K,knockouts,nMidP, numKO, timeLimit);
    end
end
warning on

%% 3. Plot envelopes
if drawEnvelope
    if ~isempty(knockouts)
        if prodMolB
            p3 = addEnv(model1,biomass, desiredProduct, knockouts, 'r', prodMol, subUptake, molarSum);
        else
            p3 = addEnv(model1,biomass, desiredProduct, knockouts, 'r');
        end

        if midPointsOn
            p={};
            for i=1:numel(ratio)-1
                color=rand(1,3);
                legendInfo{i}=['optEnvelope - Ratio= ',num2str(ratio(i))];
                if prodMolB
                    p{i} = addEnv(model1, biomass, desiredProduct, midknockouts{i}, color, prodMol, subUptake, molarSum);
                else
                    p{i} = addEnv(model1, biomass, desiredProduct, midknockouts{i}, color);
                end
            end
            try
                legend([[p1 p3],p{:}],[{'Wild-type','optEnvelope - Primary Envelope'},legendInfo])
            catch
                legend([p1 p3],{'Wild-type','optEnvelope - Primary Envelope'})
            end
        else
            legend([p1 p3],{'Wild-type','optEnvelope'})
            midknockouts = [];
        end
    else
        disp('No envelope found')
    end
    hold off
end

for i=1:length(knockouts)
    if contains(knockouts{i},'_f') || contains(knockouts{i},'_r') || contains(knockouts{i},'_b')
    knockouts{i} = knockouts{i}(1:end-2);
    end
end
