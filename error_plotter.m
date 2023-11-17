function error_plotter(filename)

load(filename)

% semilogy((s_list + r_list)*k,optimal*ones(1,length(s_list)),'k','LineWidth',5)
% hold on
% semilogy((s_list + r_list)*k,krylov_aware_untruncated,'b--*','LineWidth',3)
% semilogy((s_list + r_list)*k,krylov_aware_truncated,'b-*','LineWidth',3)
% semilogy(2*s_list*k,randSVD_truncated,'r-*','LineWidth',3)
% semilogy(2*s_list*k,exactrandSVD_truncated*ones(1,length(s_list)),'r--*','LineWidth',3)
% legend({'Optimal','Krylov aware (untruncated)','Krylov aware (truncated)',...
%     'randSVD','exact randSVD'},'location','best')
% xlabel('$s=r$','Interpreter','latex')
% ylabel('Relative Frobenius norm error','Interpreter','latex')
% title_text = append('$k = $',num2str(k));
% title(title_text,'interpreter','latex')
% set(gca,'Fontsize',14)
% hold off

semilogy(s_list,optimal*ones(1,length(s_list)),'k','LineWidth',5)
hold on
semilogy(s_list,krylov_aware_untruncated,'b--*','LineWidth',3)
semilogy(s_list,krylov_aware_truncated,'b-*','LineWidth',3)
semilogy(s_list,randSVD_truncated,'r-*','LineWidth',3)
semilogy(s_list,exactrandSVD_truncated*ones(1,length(s_list)),'r--*','LineWidth',3)
legend({'Optimal','Krylov aware (untruncated)','Krylov aware (truncated)',...
    'randSVD','exact randSVD'},'location','best')
xlabel('$s=r$','Interpreter','latex')
ylabel('Relative Frobenius norm error','Interpreter','latex')
title_text = append('$k = $',num2str(k));
title(title_text,'interpreter','latex')
set(gca,'Fontsize',16)
hold off

print(filename,'-depsc')


end