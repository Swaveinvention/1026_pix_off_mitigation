function res = dh_pattern(H, MASK, p_period, p_inactive)
    res = H;
    for r0 = 1:p_period:size(H, 1)
        for r1 = p_period-p_inactive:p_period-1
            idx = r0 + r1;
            res(idx, :) = MASK(idx, :);
        end
    end

    for c0 = 1:p_period:size(H, 2)
        for c1 = p_period-p_inactive:p_period-1
            idx = c0 + c1;
            res(:, idx) = MASK(:, idx);
        end
    end
end