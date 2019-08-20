function [s] = data_to_storage_lite(s, cnt, mou_idx, avul, rad)
    s.mou_idx(cnt.delap) = mou_idx;
    s.avul(cnt.delap) = avul;
    s.rad(cnt.delap) = rad;
end