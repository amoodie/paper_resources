function [N, dx, T, T_int] = SimpsonInt(N, dx, T, T_int)

    T_int = dx ./ 3d0 .* (  T(1) + T(N) + 2d0.*sum(T(3:2:N-2)) + 4d0.*sum(T(2:2:N-1))  );

end