function [double_diff] = double_differences(single_diff, obs_count, nsat)

design_mat = design_matrix(obs_count, nsat);
double_diff = design_mat*single_diff;

double_diff(double_diff<-1000)=0; %nan ??

end

