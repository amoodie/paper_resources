function [dtsec] = calculate_new_dtsec(dqs, Be, Bc, dx, CFLc, con, dtsecint)
    if length(dqs) > 1
        % if the passed is a vector calculate its mean
        dqs = mean(dqs);
        Be = mean(Be);
        Bc = mean(Bc);
    end
    dQs = dqs * Bc;
    dtsec = CFLc * ( (1-con.phi) * Be * (dx / -dQs) );
    
    % place limits on dt
    if dtsec < dtsecint(1)
        dtsec = dtsecint(1);
    end
    if dtsec > dtsecint(2)
        dtsec = dtsecint(2);
    end
end