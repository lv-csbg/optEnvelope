function [ActiveRxns, var, MILP]=minActiveRxns(model,matchRev,K, var, MILP, timeLimit)
% minActiveRxns determines the minimum reactions needed to be active at
% a specific point
%INPUT
% model         COBRA model structure 
% matchRev      Matching of forward and backward reactions of a reversible
%               reaction
% K             List of reactions that cannot be selected for knockout
% var           var variable for faster calculation
% MILP          MILP variable for faster calculation
% timeLimit     Time limit for gurobi optimization
%
%OUTPUT
% ActiveRxns    List of minimum active reactions at a specific point
%
% created by    Ehsan Motamedian        09/02/2022
% modified by   Kristaps Berzins        31/10/2022

if nargin < 6
    timeLimit = inf;
end

if isempty(MILP)

    [nMetsR,nRxnsR] = size(model.S);
    U=1001;

    % Create MILP problem for first point
    MILP = struct;
    MILP.A= model.S;
    MILP.b=model.b;
    MILP.csense(1:nMetsR,1) = 'E';
    MILP.c=zeros(nRxnsR,1);
    MILP.lb=model.lb;
    MILP.ub=model.ub;
    MILP.vartype(1:nRxnsR)='C';
    n=0;
    for i=1:nRxnsR
        if ismember(i,K)==0
            n=n+1;
            MILP.A(nMetsR+n,i)=1;
            MILP.A(nMetsR+n,nRxnsR+n)=-U;
            MILP.b(nMetsR+n,1)=0;MILP.csense(nMetsR+n,1)='L';
            MILP.c(nRxnsR+n)=1;
            MILP.lb(nRxnsR+n)=0;MILP.ub(nRxnsR+n)=1;
            MILP.vartype(nRxnsR+n)='B';
            var(n)=i;
        end
    end
    k=0;
    for or=1:length(matchRev)
        if matchRev(or)~=0
            k=k+1;
            MILP.A(nMetsR+n+k,nRxnsR+find(ismember(var,or)))=1;
            MILP.A(nMetsR+n+k,nRxnsR+find(ismember(var,matchRev(or))))=1;
            MILP.b(nMetsR+n+k,1)=1;
            MILP.csense(nMetsR+n+k,1)='L';
            matchRev(matchRev(or))=0;
        end
    end
    MILP.x0=zeros(nRxnsR+n,1);
    MILP.osense=1;

    % Solve MILP problem for first point
    try
        s=solveCobraMILP(MILP, 'timeLimit', timeLimit);
        ActiveRxns=model.rxns(var(s.int==1));
    catch
        ActiveRxns=[];
    end

else
    % Solve MILP problem for next point if previous point failed
    MILP.lb(1:length(model.lb))=model.lb;
    MILP.ub(1:length(model.ub))=model.ub;
    try
        s=solveCobraMILP(MILP, 'timeLimit', timeLimit);
        ActiveRxns=model.rxns(var(s.int==1));
    catch
        ActiveRxns=[];
    end
end


