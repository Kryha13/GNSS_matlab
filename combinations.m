function [phase, code] = combinations(comb)

switch comb                                                
    case "com1" %GPS
        L1 = "L7Q";
        L2 = "L5Q";
        P1 = "C7Q";
        P2 = "C5Q";
    case "com2" %GLONASS
        L1 = "L6C";
        L2 = "L7Q";
        P1 = "C6C";
        P2 = "C7Q";
    case "com3" %Galileo
        L1 = "L1C";
        L2 = "L6C";
        P1 = "C1C";
        P2 = "C6C";
end
    phase = [L1 L2];
    code = [P1 P2];
end
