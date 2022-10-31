clc
clear
% 
changeCobraSolver('gurobi','all',0);
model=readCbModel('iJR904');
biomass='BIOMASS_Ecoli';
desiredProduct='EX_ac_e';
desiredProductName='Acetate';
[model,matchRev] = convertToIrreversible(model);
K=[find(contains(model.rxns,'EX_')); findRxnIDs(model, 'H2Ot_f'); findRxnIDs(model, 'H2Ot_b')];
knockouts = optEnvelope(model,biomass,desiredProduct,desiredProductName,matchRev,K, 60)
