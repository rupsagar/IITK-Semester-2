function [NotConv,ChkVal0] = ConvTest(ConvCrit,tol,ChkVal0,NotConv,Index,dr,Ru)

if strcmp(ConvCrit,"Relative Energy Increment")
    if Index==0
        ChkVal0 = abs(Ru'*dr);
    else
        StrChk = abs(Ru'*dr)/ChkVal0;
        if StrChk<tol
            NotConv = false;
        end
    end
elseif strcmp(ConvCrit,"Energy Increment")
    StrChk = abs(Ru'*dr);
    if StrChk<tol
        NotConv = false;
    end
    ChkVal0 = [];
end
    
end