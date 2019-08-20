function [mou_idx] = check_for_mou_prog(eta, zed, depthick, Hform, rad0_idx, mou0_idx, mouSwitch)
    switch mouSwitch
        case 'bartip'
            if any(eta(mou0_idx:end) > zed(mou0_idx:end)-Hform)
                [~, max_idx] = max(depthick);
                adj = 2;
                if ~ (max_idx - adj <= mou0_idx)
                    mou_idx = max_idx - adj;
                else
                    mou_idx = mou0_idx;
                end
            else
                mou_idx = mou0_idx;
            end
        case 'na'
            mou_idx = mou0_idx;
        case 'endbar' % just pick the end of the bar
            if any(depthick(mou0_idx:end) > depthick(mou0_idx)) 
                mou_idx = mou0_idx + 1;
            else
                mou_idx = mou0_idx;
            end
        case 'endbar+above' % just pick the end of the bar
            above_idx = eta(rad0_idx:end) > zed(rad0_idx:end)-Hform;
            if any(above_idx)
                mou_idx = rad0_idx + find(above_idx, 1, 'last') - 1;
            else
                mou_idx = mou0_idx;
            end
        case 'endbar+dep'
            if and( any(depthick(mou0_idx:end) > depthick(mou0_idx)), any(Hbasin(mou0_idx:end) - depthick(mou0_idx:end) < Hform) )
                mou_idx = mou0_idx + 1;
            else
                mou_idx = mou0_idx;
            end
    end
end