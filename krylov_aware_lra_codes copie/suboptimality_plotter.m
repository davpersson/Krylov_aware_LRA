function suboptimality_plotter(filename)

load(filename)

semilogy(s_list,krylov_aware_,'b-*','LineWidth',3)
hold on
semilogy(s_list,subspace_iteration,'r-*','LineWidth',3)

legend({'Krylov aware','Subspace iteration with Lanczos'},'location','best')
xlabel('$s=r$','Interpreter','latex')
ylabel('$\varepsilon$','Interpreter','latex')
title_text = append('$k = $',num2str(k));
title(title_text,'interpreter','latex')
set(gca,'Fontsize',16)
hold off
hold on
for i = 1:length(s_list)
    text(s_list(i),2*subspace_iteration(i),append('q=',num2str(q_optimal(i))),'Color','k','FontSize',12)
end
hold off

print(filename,'-depsc')


end