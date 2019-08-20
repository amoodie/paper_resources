function [s] = data_to_storage(s, cnt, eta, H, zed, mou_idx, avul, rad, Qw, Lb, Be, Bc, qs)
    s.eta(:, cnt.delap) = eta;
    s.H(:, cnt.delap) = H;
    s.zed(:, cnt.delap) = zed;
    s.mou_idx(cnt.delap) = mou_idx;
    s.avul(cnt.delap) = avul;
    s.rad(cnt.delap) = rad;
    s.Qw(cnt.delap) = Qw;
    s.Lb(cnt.delap) = Lb;
    s.Be(:, cnt.delap) = uint32(Be);
    s.Bc(:, cnt.delap) = uint32(Bc);
    s.qs(:, cnt.delap) = qs;
end