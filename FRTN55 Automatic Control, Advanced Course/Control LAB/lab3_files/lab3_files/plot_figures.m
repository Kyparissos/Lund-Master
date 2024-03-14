line_styles = {'--b','-b','--r','-r'}
figure(1)
if length(tout) == length(Left_figure)
    clf 
    hold off
    hold on
    for k = 1:4
        plot(tout, Left_figure(:,k), line_styles{k},...
            'LineWidth',2)
    end
    title('Left Figure')
    ylim([0 16])
    legend('r_1', 'h_1', 'r_2', 'h_2')
end
figure(2)
if length(tout) == length(Right_figure)
    clf 
    hold off
    hold on
    for k = 1:4
        plot(tout, Right_figure(:,k), line_styles{k},...
            'LineWidth',2)
    end
    title('Right Figure')
    ylim([0 16])
    legend('r_1', 'h_1', 'r_2', 'h_2')
end
