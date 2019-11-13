function [C,damp_ratios] = damp_mat(K_S,M_S,w_tr)

A = [
    .5*w_tr(1), .5/w_tr(1);
    .5*w_tr(2), .5/w_tr(2);
    ];

B = [ .01; .01];

sol = A\B;
C = sol(1) * K_S + sol(2) * M_S;


damp_ratios = zeros(1,length(w_tr));

fprintf('Damping ratios under proportional assumption : \n');
for i = 1:length(w_tr)
   ratio = .5*(sol(1)*w_tr(i)+sol(2)/w_tr(i)); 
   damp_ratios(i) = ratio;
   fprintf('Mode %i : \x03B5 = %.2f %%\n',i,100*ratio);
end
end