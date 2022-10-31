function [ActiveRxns,var,s]=minActiveRxns(model,matchRev,K)
% minActiveRxns determines the minimum reactions needed to be active at
% a specific point
%INPUT
% model         COBRA model structure 
% matchRev      Matching of forward and backward reactions of a reversible
%               reaction
% K             List of reactions that cannot be selected for knockout
%
%OUTPUT
% ActiveRxns     List of minimum active reactions at a specific point
%
% created by    Ehsan Motamedian        09/02/2022
% modified by   Kristaps Berzins        31/10/2022

[nMetsR,nRxnsR] = size(model.S);
U=1001;

% Create MILP problem
model2.A= model.S;
model2.b=model.b;
model2.csense(1:nMetsR,1) = 'E';
model2.c=zeros(nRxnsR,1);
model2.lb=model.lb;
model2.ub=model.ub;
model2.vartype(1:nRxnsR)='C';
n=0;
for i=1:nRxnsR
    if ismember(i,K)==0
        n=n+1;
        model2.A(nMetsR+n,i)=1;
        model2.A(nMetsR+n,nRxnsR+n)=-U;
        model2.b(nMetsR+n,1)=0;model2.csense(nMetsR+n,1)='L';
        model2.c(nRxnsR+n)=1;
        model2.lb(nRxnsR+n)=0;model2.ub(nRxnsR+n)=1;
        model2.vartype(nRxnsR+n)='B';
        var(n)=i;
    end
end
k=0;
for or=1:length(matchRev)
    if matchRev(or)~=0
        k=k+1;
        model2.A(nMetsR+n+k,nRxnsR+find(ismember(var,or)))=1;
        model2.A(nMetsR+n+k,nRxnsR+find(ismember(var,matchRev(or))))=1;
        model2.b(nMetsR+n+k,1)=1;
        model2.csense(nMetsR+n+k,1)='L';
        matchRev(matchRev(or))=0;
    end
end
model2.x0=zeros(nRxnsR+n,1);
model2.osense=1;

% Solve MILP problem
try
    s=solveCobraMILP(model2);
    ActiveRxns=model.rxns(var(s.int==1));
catch
    ActiveRxns=[];
end


