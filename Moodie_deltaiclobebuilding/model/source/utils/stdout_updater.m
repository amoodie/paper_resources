function stdout_updater(cnt)
    % could add support to output more information in the future...

    disp(['n/N: ' num2str(cnt.n) '/' num2str(cnt.N) '; Elapsed time: ' num2str(cnt.y) 'y ' num2str(cnt.d-1) 'd'])
    % disp(['Qw: ', num2str(Qw(cnt.delap+1)), ', dt: ', num2str(dtsec)])
    % disp(['dt = ' num2str(dtsec / 86400) ' of a day'])
    % disp(['s_max = 10^', num2str(log10(max(S)))])

end

