function [B] = B_matrix(obs_number, phase_freqencies)
             
B0 = zeros(length(phase_freqencies)*obs_number, length(phase_freqencies)*obs_number);
B00 = zeros(obs_number, obs_number);
for i=1:length(phase_freqencies)
    
   B1{i} = eye(obs_number)*phase_freqencies(i);
    
end

B = [B1(1) B00 B00 B00 B00;
    B00 B1(2) B00 B00 B00;
    B00 B00 B1(3) B00 B00;
    B00 B00 B00 B1(4) B00;
    B00 B00 B00 B00 B1(5);
    ];

B = cell2mat(B);
B = [B;B0];

% B = [];
% for ii=1:length(B1)
%   B=blkdiag(B, B1(1,ii));
%   
% end
% B11 = blkdiag(B1(:));

% B = [B11; B0];

end

