function svk_error_plotter(filename)

% Plots results from svk_error_test.m

load(filename)

semilogy(2*k*s_list,optimal*ones(1,length(s_list)),'k','LineWidth',5)
hold on
semilogy(2*k*s_list,krylov_aware_untruncated,'b--*','LineWidth',3,'MarkerSize', 12)
semilogy(2*k*s_list,krylov_aware_truncated,'b-+','LineWidth',3,'MarkerSize', 15)
semilogy(2*k*s_list,svk_krylov_aware_untruncated,'m-|','LineWidth',3,'MarkerSize', 12)
semilogy(2*k*s_list,svk_krylov_aware_truncated,'m--^','LineWidth',3,'MarkerSize', 15)
legend({'Optimal','Krylov aware (untruncated)','Krylov aware (truncated)',...
    'SV Krylov aware (untruncated)','SV Krylov aware (truncated)'},'location','northeast','interpreter','latex')
xlabel('Matrix vector products with $\mathbf{A}$','Interpreter','latex')
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