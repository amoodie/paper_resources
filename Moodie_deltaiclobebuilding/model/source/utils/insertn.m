function y = insertn(p, q, n)
    % adapted from code at
    % http://stackoverflow.com/questions/14631182/padding-an-array-after-every-nth-element
    % insert q into p at every n
    N = numel(p);
    p_pad = [p(:); zeros((n - mod(N, n)) * (mod(N, n) > 0), 1)];
    y = [reshape(p_pad, n, []); repmat(q(:), 1, numel(p_pad) / n)];
    y = y(1:N + numel(q) * fix(N / n));
    if isrow(p)
        y = y';
    end
end