function [FIG_r0] = plotR0Data(data15, data16)
    col15 = [0 0 1];
    col16 = [1 0 0];
    FIG_r0 = figure();
    subplot(2, 3, 1)
        hold on
        plot(data15.modelEvalCsN', data15.modelEvalZsN', 'Color', col15)
        plot(data16.modelEvalCsN', data16.modelEvalZsN', 'Color', col16)
        text(0.05, 0.2, 'Zr = ws/\kappa u_*')
        text(0.05, 0.1, 'solve for u_*')
        xlabel('conc relative (Rousean coords)')
        ylabel('z relative (Rousean coords)')
    subplot(2, 3, 2)
        cla
        hold on
        plot(data15.lotw.Us', data15.modelEvalZs', 'Color', col15)
        plot(data16.lotw.Us', data16.modelEvalZs', 'Color', col16)
        text(2, 4.5, 'U_z = (u_* / 0.41) * log(z/z0)')
        text(2, 3.5, '$C_{HU} = \int_0^H (u_z c_z dz)/(H\bar{u})$', 'Interpreter', 'latex')
        xlabel('velocity U (m/s)')
        ylabel('elev z (m)')
    subplot(2, 3, 3)
        hold on
        plot(data15.C_HU, data15.fit.cb, 'Marker', 'o', 'LineStyle', 'none', 'Color', col15)
        plot(data16.C_HU, data16.fit.cb, 'Marker', 'o', 'LineStyle', 'none', 'Color', col16)
        xlabel('C_{HU}')
        ylabel('cb')
    subplot(2, 3, 4)
        hold on
        plot(data15.fit.ustar, data15.r0, 'Marker', 'o', 'LineStyle', 'none', 'Color', col15)
        plot(data16.fit.ustar, data16.r0, 'Marker', 'o', 'LineStyle', 'none', 'Color', col16)
        xlabel('u_*')
        ylabel('r0')
    subplot(2, 3, 5)
        loglog(data15.Zu, data15.fit.cb, 'Marker', 'o', 'LineStyle', 'none', 'Color', col15)
        hold on
        plot(data16.Zu, data16.fit.cb, 'Marker', 'o', 'LineStyle', 'none', 'Color', col16)
        xlabel('Zu')
        ylabel('cb')
    subplot(2, 3, 6)
        plot(data15.fit.ustar, data15.fit.cb, 'Marker', 'o', 'LineStyle', 'none', 'Color', col15);
        hold on
        plot(data16.fit.ustar, data16.fit.cb, 'Marker', 'o', 'LineStyle', 'none', 'Color', col16)
        xlabel('u_*')
        ylabel('cb')
end