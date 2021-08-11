function [wt,x] = ParamIntgr(IntgrType,MemL,NumMem,NumIntgrPts)

xi = cell(NumMem,1);
wt = cell(NumMem,1);
x = cell(NumMem,1);

for mem = 1:NumMem
    if strcmp(IntgrType(mem),"Legendre")    
        switch NumIntgrPts(mem)
            case 1
                xi{mem} = 0;
                wt{mem} = 2;
            case 2
                xi{mem} = 0.577350269*[-1,1];
                wt{mem} = [1,1];
            case 3
                xi{mem} = 0.77459669*[-1,0,1];
                wt{mem} = [5/9,8/9,5/9];
            case 4
                xi{mem} = [-0.861136312,-0.339981044,0.339981044,0.861136312];
                wt{mem} = [0.3478548,0.6521452,0.6521452,0.3478548];
            case 5
                xi{mem} = [-0.906179846,-0.538469310,0,0.538469310,0.906179846];
                wt{mem} = [0.2369269,0.538469310,0.5688889,0.538469310,0.2369269];
        end
    elseif strcmp(IntgrType(mem),"Lobatto")
        switch NumIntgrPts(mem)
            case 2
                xi{mem} = [-1,1];
                wt{mem} = [1,1];
            case 3
                xi{mem} = [-1,0,1];
                wt{mem} = [1/3,4/3,1/3];
            case 4
                xi{mem} = [-1,-1/sqrt(5),1/sqrt(5),1];
                wt{mem} = [1/6,5/6,5/6,1/6];
            case 5
                xi{mem} = [-1,-sqrt(3/7),0,sqrt(3/7),1];
                wt{mem} = [0.1,49/90,32/45,49/90,0.1];
            case 6
                xi{mem} = [-1,-sqrt(1/3+2*sqrt(7)/21),-sqrt(1/3-2*sqrt(7)/21),...
                    sqrt(1/3-2*sqrt(7)/21),sqrt(1/3+2*sqrt(7)/21),1];
                wt{mem} = [1/15,(14-sqrt(7))/30,(14+sqrt(7))/30,(14+sqrt(7))/30,(14-sqrt(7))/30,1/15];
            case 7
                xi{mem} = [-1,-sqrt(5/11+2/11*sqrt(5/3)),-sqrt(5/11-2/11*sqrt(5/3)),0,...
                    sqrt(5/11-2/11*sqrt(5/3)),sqrt(5/11+2/11*sqrt(5/3)),1];
                wt{mem} = [1/21,(124-7*sqrt(15))/350,(124+7*sqrt(15))/350,256/525,...
                    (124+7*sqrt(15))/350,(124-7*sqrt(15))/350,1/21];                
        end
    end
    x{mem} = MemL(mem)/2*(1+xi{mem});
end

end