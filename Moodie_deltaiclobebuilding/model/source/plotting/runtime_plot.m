function [fig] = runtime_plot(fig, s, x, eta, zed, H, Hnform, U, Qw, mou_idx, rad_idx, Bc, Be, dx, cnt)
    figure(fig)
    subplot(5, 1, 1:2)
        cla; hold on;
        plot(x', s.eta(:,1), '-', 'LineWidth', 1.5, 'Color', [0.5 0.5 0.5])
        plot(x, eta, 'r-', 'LineWidth', 1.5)
        plot(x, zed - Hnform, 'k:', 'LineWidth', 1.5)
        plot(x, zed, 'g-', 'LineWidth', 1.5)
        plot(x, eta+H, 'b-', 'LineWidth', 1.5)
        plot([s.mou_idx(1)*dx s.mou_idx(1)*dx], ylim, 'k--', 'LineWidth', 1.5)
        plot([mou_idx*dx mou_idx*dx], ylim, 'k-', 'LineWidth', 1.5)
    subplot(5, 1, 3)
        cla; hold on;
        plot(x, U, 'LineWidth', 1.5)
        plot([s.mou_idx(1)*dx s.mou_idx(1)*dx], ylim, 'k--', 'LineWidth', 1.5)
        plot([mou_idx*dx mou_idx*dx], ylim, 'k-', 'LineWidth', 1.5)
        xlim([0, max(x)])
        ylabel('U')
    subplot(5, 1, 4)
        cla; hold on;
        plot(x, Bc, 'LineWidth', 1.5)
        plot(x, Be, 'LineWidth', 1.5)
        plot([s.mou_idx(1)*dx s.mou_idx(1)*dx], [0, max(Be)], 'k--', 'LineWidth', 1.5)
        plot([mou_idx*dx mou_idx*dx], [0, max(Be)], 'k-', 'LineWidth', 1.5)
        ylabel('width')
    subplot(5, 1, 5)
        cla; hold on;
        plot(1:cnt.delap+1, Qw(1:cnt.delap+1), 'LineWidth', 1.5)
    drawnow
end

