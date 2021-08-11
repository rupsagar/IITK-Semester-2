function [NumTotalDataPts,kIni,e,rx,Ssi,vsi,Ssk,vsk,Kv,Smk,vmk,Qk,qk,Ks,Rk,rk,i] = ...
    Initialize(NumSecDOF,NumSec,NumFib,NumBasDOF,NumStrDOF,NumMem,NumMemDOF,Kv0,Ks0,State,NumDataPts)

NumDataPtsPrev = State.NumTotalDataPts;

if NumDataPtsPrev==0
    NumTotalDataPts = NumDataPts;
    kIni = 2;
elseif NumDataPtsPrev~=0
    NumTotalDataPts = NumDataPtsPrev+NumDataPts-1;
    kIni = NumDataPtsPrev+1;
end

e = cell(NumMem,1);
rx = cell(NumMem,1);
Ssi = cell(NumMem,1);
vsi = cell(NumMem,1);
Ssk = cell(NumMem,1);
vsk = cell(NumMem,1);

Kv = zeros(NumBasDOF,NumBasDOF,NumMem,NumTotalDataPts);
Smk = zeros(NumBasDOF,NumMem,NumTotalDataPts);
vmk = zeros(NumBasDOF,NumMem,NumTotalDataPts);
Qk = zeros(NumMemDOF,NumMem,NumTotalDataPts);
qk = zeros(NumMemDOF,NumMem,NumTotalDataPts);
Ks = zeros(NumStrDOF,NumStrDOF,NumTotalDataPts);
Rk = zeros(NumStrDOF,NumTotalDataPts);
rk = zeros(NumStrDOF,NumTotalDataPts);

i = zeros(NumTotalDataPts,1);

% INITIALIZATION FOR 1ST ITERATION
for mem = 1:NumMem
    e{mem} = zeros(NumFib,NumSec(mem));
    rx{mem} = zeros(NumSecDOF,NumSec(mem));   
    Ssi{mem} = zeros(NumSecDOF,NumSec(mem));
    vsi{mem} = zeros(NumSecDOF,NumSec(mem));
    Ssk{mem} = zeros(NumSecDOF,NumSec(mem),NumTotalDataPts);
    vsk{mem} = zeros(NumSecDOF,NumSec(mem),NumTotalDataPts);
end

if NumDataPtsPrev==0
    for mem = 1:NumMem
        Kv(:,:,mem,1) = Kv0(:,:,mem);
    end
    Ks(:,:,1) = Ks0;
elseif NumDataPtsPrev>0
    for mem = 1:NumMem
        Ssk{mem}(:,:,1:NumDataPtsPrev) = State.Ssk{mem}(:,:,:);
        vsk{mem}(:,:,1:NumDataPtsPrev) = State.vsk{mem}(:,:,:);
    end
    Kv(:,:,:,1:NumDataPtsPrev) = State.Kv;
    Smk(:,:,1:NumDataPtsPrev) = State.Smk;
    vmk(:,:,1:NumDataPtsPrev) = State.vmk;
    Qk(:,:,1:NumDataPtsPrev) = State.Qk;
    qk(:,:,1:NumDataPtsPrev) = State.qk;
    Ks(:,:,1:NumDataPtsPrev) = State.Ks;
    Rk(:,1:NumDataPtsPrev) = State.Rk;
    rk(:,1:NumDataPtsPrev) = State.rk;
end

end