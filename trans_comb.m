function [phase, code] = trans_comb(comb)

switch comb                                                
    case "trans1" %BeiDou
        L1 = "L6C";
        L2 = "L7Q";
        L3 = "L5Q";
        P1 = "C6C";
        P2 = "C7Q";
        P3 = "C5Q";
    case "trans2" %QZSS PRN-192
        L1 = "L1C";
        L2 = "L7Q";
        L3 = "L5Q";
        P1 = "C1C";
        P2 = "C7Q";
        P3 = "C5Q";
end
    phase = [L1 L2 L3];
    code = [P1 P2 P3];
end