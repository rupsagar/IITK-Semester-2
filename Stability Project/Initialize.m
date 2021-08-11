function Init = Initialize(NumNodeDOF,ElementType,XY,MemCon,StrBC,LoadPattern,NumStep,E,A,I)

% BASIC PARAMETERS
Init.NumNode = numel(XY)/2;
Init.NumMem = size(MemCon,1);
Init.NumMemDOF = 2*NumNodeDOF;
Init.NumTotalDOF = NumNodeDOF*Init.NumNode;
Init.NumStrDOF = sum(~StrBC(:));


% INITIALIZATION
Init.MemDOF = zeros(Init.NumMem,Init.NumMemDOF);
Init.L0 = zeros(Init.NumMem,1);
Init.LoadFac = zeros(NumStep,1);
Init.Ks = zeros(Init.NumStrDOF,Init.NumStrDOF,NumStep);
Init.RefLoad = zeros(Init.NumTotalDOF,1);
Init.Re = zeros(Init.NumStrDOF,NumStep);
Init.Ru = zeros(Init.NumStrDOF,1);
Init.drSys = zeros(Init.NumTotalDOF,1);
Init.XYk = zeros(2*Init.NumNode,NumStep);
Init.ii = zeros(NumStep,1);
Init.Beta0 = zeros(Init.NumMem,1); % REQUIRED FOR FRAME ELEMENTS
Init.Thetak = zeros(Init.NumNode,NumStep); % REQUIRED FOR FRAME ELEMENTS

Init.BCMap = ~StrBC(:);

% Init.KBCMap = logical(Init.BCMap'.*Init.BCMap);
Init.KBCMap = false(numel(Init.BCMap));
for iii = 1:numel(Init.BCMap)
    for jjj = 1:numel(Init.BCMap)
        Init.KBCMap(iii,jjj) = Init.BCMap(iii)*Init.BCMap(jjj);
    end
end

Init.StrXYTheta = [false(2,Init.NumNode);true(1,Init.NumNode)];
Init.XYThetaMap = ~Init.StrXYTheta(:);

NumLoadPattern = size(LoadPattern,1);
for i = 1:NumLoadPattern
    LoadNodeID = LoadPattern(i,1);
    Init.RefLoad(NumNodeDOF*LoadNodeID-(NumNodeDOF-1):NumNodeDOF*LoadNodeID) = LoadPattern(i,2:NumNodeDOF+1)';
end
Init.RefLoad = Init.RefLoad(Init.BCMap);

for mem = 1:Init.NumMem
    MemNode = MemCon(mem,:);
%     if NumNodeDOF==2
%         Init.MemDOF(mem,:) = [2*MemNode(1)-1,2*MemNode(1),2*MemNode(2)-1,2*MemNode(2)];
%     elseif NumNodeDOF==3
        Init.MemDOF(mem,:) = [3*MemNode(1)-2,3*MemNode(1)-1,3*MemNode(1),...
            3*MemNode(2)-2,3*MemNode(2)-1,3*MemNode(2)];
%     end
    Node1XY = XY([2*MemNode(1)-1,2*MemNode(1)]);
    Node2XY = XY([2*MemNode(2)-1,2*MemNode(2)]);
    ri = Node2XY-Node1XY;
    Init.L0(mem) = norm(ri);
    Init.Beta0(mem) = atan(ri(2)/ri(1));
end

Init.XYk(:,1) = XY;
[Init.Ks(:,:,1),~] = CorAssemble(ElementType,MemCon,Init.XYk(:,1),Init.Thetak(:,1),Init.L0,Init.Beta0,...
    Init.NumMem,Init.NumTotalDOF,Init.NumStrDOF,Init.NumMemDOF,Init.MemDOF,Init.KBCMap,Init.BCMap,E,A,I);