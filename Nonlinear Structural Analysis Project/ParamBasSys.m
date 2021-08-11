function [bx,Kv0,avq] = ParamBasSys(BasicSystem,NumBasDOF,NumMem,MemL,ASec,IzSec,MatModel,apq)

if strcmp(MatModel{1},"Bilinear")
    E = MatModel{3};
elseif strcmp(MatModel{1},"GMP")
    E = MatModel{3};
elseif strcmp(MatModel{1},"Elastic")
    E = MatModel{2};
end

bx = cell(NumMem,1);
avp = zeros(NumBasDOF,6,NumMem);
Kv0 = zeros(NumBasDOF,NumBasDOF,NumMem);
avq = zeros(NumBasDOF,6,NumMem);

if strcmp(BasicSystem,"Cantilever")
    for mem = 1:NumMem
        bx{mem} = @(x)[1,0,0;0,MemL(mem)-x,1];
        avp(:,:,mem) = [-1,0,0,1,0,0;
                      0,-1,-MemL(mem),0,1,0;
                      0,0,-1,0,0,1];        
        Kv0(:,:,mem) = [ASec{mem}(1)*E/MemL(mem),0,0;
                        0,12*E*IzSec(mem)/MemL(mem)^3,-6*E*IzSec(mem)/MemL(mem)^2;
                        0,-6*E*IzSec(mem)/MemL(mem)^2,4*E*IzSec(mem)/MemL(mem)];
        avq(:,:,mem) = avp(:,:,mem)*apq(:,:,mem);
    end
elseif strcmp(BasicSystem,"Simply Supported")||strcmp(BasicSystem,"Simply supported")
    for mem = 1:NumMem
        bx{mem} = @(x)[1,0,0;0,(x/MemL(mem)-1),x/MemL(mem)];
        avp(:,:,mem) = [-1,0,0,1,0,0;
                      0,1/MemL(mem),1,0,-1/MemL(mem),0;
                      0,1/MemL(mem),0,0,-1/MemL(mem),1];
        Kv0(:,:,mem) = [ASec{mem}(1)*E/MemL(mem),0,0;
                        0,4*E*IzSec(mem)/MemL(mem),2*E*IzSec(mem)/MemL(mem);
                        0,2*E*IzSec(mem)/MemL(mem),4*E*IzSec(mem)/MemL(mem)];
        avq(:,:,mem) = avp(:,:,mem)*apq(:,:,mem);
    end      
end

end