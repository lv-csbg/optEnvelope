clc
clear
% Choose solver
changeCobraSolver('gurobi','all',0);

% Load model
model=readCbModel('iJR904');

% Describe biomass
biomass='BIOMASS_Ecoli';

% Describe target product
desiredProduct='EX_ac_e';
desiredProductName='Acetate';
prMolarMass = 60;

% Convet model to irreversible form
[model,matchRev] = convertToIrreversible(model);
% Create a list of protected reactions (their IDs)
K=[find(contains(model.rxns,'EX_')); findRxnIDs(model, 'H2Ot_f'); findRxnIDs(model, 'H2Ot_b')];

%run optEnvelope
knockouts = optEnvelope(model,biomass,desiredProduct,desiredProductName,matchRev,K, prMolarMass)
