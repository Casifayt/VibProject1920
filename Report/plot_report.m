% t_res = [197.258,294.693,346.871,413.507,428.49,1050.147,2101.075,3508.331,7790.724];
% t_discr = [24.289,42.597,53.757,63.702,54.511,73.359,78.455,97.113,119.379];
% freq_newmark_conv = [0.02, 0.2926829268,0.399024938,0.399900025,0.399900002,0.399999];
% N_newmrk_conv = [10,1,.1,.01,.001,.0001];
% freq1_newmark_conv = [2.1994501375, 2.1999450014, 2.1999945];
% N1_newmark_conv = [.01,.001,.0001];
% freq2_newmark_conv = [1.2219451372, 1.2746813297, 1.2749681258, 1.2749968125];
% N2_newmark_conv = [.1,.01,.001,.0001];
% t_comp_newmark_conv = [0.11,0.17,0.81,6.91,68.38,798.51];
% tcomp_CB = [0.10961,0.140497,0.15954,0.165064,0.176534,0.203152,0.231173,0.371916];

e_rel_CB_24 = [120.35,16.13,.55,.14,.14,.11,.037,.03,.017,.01,.0089,.0083,.0059];
e_rel_CB_48 = [49.78,1.44,.43,.28,.28,.21,.16,0.03,.021,.018,.018,.015,.014];
e_rel_CB_72 = [14.93,1.44,.4,.4,.36,.3,.12,.03,.029,.027,.026,.022,.021];
N_e = [1,5,10,15,20,30,40,50,60,70,80,90,100];
Nt = [1,5,10,15,20,30,50,100];

fig = figure;
set(fig,'defaultAxesColorOrder',[[0.8500 0.3250 0.0980] ;[0 0.4470 0.7410] ]);

yyaxis left
plot(N_e,e_rel_CB_24, '.-', 'Color',[0.9 0.4 0.2]); hold on;
plot(N_e,e_rel_CB_48, '-', 'Color','r');
plot(N_e,e_rel_CB_72, '--', 'Color',[0.6350 0.0780 0.1840]);
% set(gca,'xscale','log');  
% set(gca,'XDir','Reverse');
set(gca,'yscale','log');  

ylabel('\epsilon_{rel}^{avg} [%]'); xlabel('Number of modes');



% 
yyaxis right;
plot(Nt,tcomp_CB,'x', 'MarkerSize',8);  grid on; 
set(gca,'yscale','log');
ylabel('Time [s]');
legend('24 DOFs', '48 DOFs', '72 DOFs','Averaged computation time');



% xlabel('Nbr of element');
% legend({'Time of resolution', 'Average relative error'},'Location','north');
