# optEnvelope
## 1. Initialize cobra toolbox and set solver
initCobraToolbox(false) <br />
changeCobraSolver('gurobi','all',0)

## 2. Load model and set desired product
load('iJR904.mat') <br />
desiredProduct='EX_ac_e';

## 3. Run optEnvelope
[main] = optEnvelope(iJR904, desiredProduct,'timeLimit', 600);

##3.1 Run optEnvelope with mid envelopes
[main, mid] = optEnvelope(iJR904, desiredProduct, 'midPoints', 10)

## 4. Preprint
https://www.biorxiv.org/content/10.1101/2023.03.10.532079v1
