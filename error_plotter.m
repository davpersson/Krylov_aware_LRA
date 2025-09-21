function error_plotter(filename)

% Plots the results from error_test.m

load(filename)

semilogy(s_list,optimal*ones(1,length(s_list)),'k','LineWidth',5)
hold on
semilogy(s_list,krylov_aware_untruncated,'b--*','LineWidth',3,'MarkerSize', 12)
semilogy(s_list,krylov_aware_truncated,'b-+','LineWidth',3,'MarkerSize', 12)
semilogy(s_list,randSVD_truncated,'r-o','LineWidth',3,'MarkerSize', 12)
semilogy(s_list,exactrandSVD_truncated*ones(1,length(s_list)),'r--s','LineWidth',3,'MarkerSize', 12)
legend({'Optimal','Krylov aware (untruncated)','Krylov aware (truncated)',...
    'randSVD','exact randSVD'},'location','best','Interpreter','latex')
xlabel('$s=r$','Interpreter','latex')
ylabel('Relative Frobenius norm error','Interpreter','latex')
title_text = append('$k = $',num2str(k));
title(title_text,'interpreter','latex')
curtick = get(gca, 'xTick');
xticks(unique(round(curtick)));
set(gca,'Fontsize',18)
set(gca,'TickLabelInterpreter','latex');
hold off

print(filename,'-dpng')


end