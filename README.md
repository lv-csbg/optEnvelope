# optEnvelope
## 1. Initialize cobra toolbox and set solver
initCobraToolbox(false)
changeCobraSolver('gurobi','all',0)

## 2. Load model and set desired product
load('iJR904.mat')
desiredProduct='EX_ac_e';

## 3. Run optEnvelope
[knockouts, midknockouts] = optEnvelope(iJR904, desiredProduct,'timeLimit', 600);
