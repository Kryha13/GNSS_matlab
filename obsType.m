function [obs_type1, obs_type2] = obsType(obs_name)

C1a = []; 
C1b = [];
L1a = [];
L1b = [];
C6a = [];
C6b = [];
L6a = [];
C5a = [];
C5b = [];
L5a = [];
L5b = [];
C7a = [];
C7b = [];
L7a = [];
L7b = [];
C8a = [];
C8b = [];
L8a = [];
L8b = [];


switch obs_name                                                
    case 'C1C'
        obs_type1 = C1a;
        obs_type2 = C1b;
    case 'L1C'
        obs_type1 = L1a;
        obs_type2 = L1b;
    case 'C6C' 
        obs_type1 = C6a;
        obs_type2 = C6b;
    case 'L6C' 
        obs_type1 = L6a;
        obs_type2 = L6b;
    case 'C5Q'
        obs_type1 = C5a;
        obs_type2 = C5b;
    case 'L5Q' 
        obs_type1 = L5a;
        obs_type2 = L5b;
    case 'C7Q' 
        obs_type1 = C7a;
        obs_type2 = C7b;
    case 'L7Q' 
        obs_type1 = L7a;
        obs_type2 = L7b;
    case 'C8Q' 
        obs_type1 = C8a;
        obs_type2 = C8b;
    case 'L8Q'
        obs_type1 = L8a;
        obs_type2 = L8b;
end

end