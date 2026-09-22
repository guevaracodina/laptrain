function y = replaceNanBlocks(y)
    assert(isvector(y), 'Input must be a vector');
    wasRow = isrow(y);
    y = y(:);
    n = numel(y);

    isn_vec = isnan(y);
    good = ~isn_vec;
    firstGood = find(good,1,'first');
    lastGood  = find(good,1,'last');
    if isempty(firstGood)
        error('All values are NaN; cannot fill.');
    end
    % handle leading NaNs
    if firstGood > 1
        y(1:firstGood-1) = y(firstGood);
    end
    % handle trailing NaNs
    if lastGood < n
        y(lastGood+1:end) = y(lastGood);
    end

    % recompute NaN mask
    isn_vec = isnan(y);
    d = diff([0; isn_vec; 0]);
    run_starts = find(d == 1);
    run_ends   = find(d == -1) - 1;

    for k = 1:numel(run_starts)
        i0 = run_starts(k);
        i1 = run_ends(k);
        len = i1 - i0 + 1;
        v_before = y(i0 - 1);
        v_after  = y(i1 + 1);

        if len == 1
            % simple average if only one NaN
            y(i0) = (v_before + v_after)/2;
        else
            % make linear interpolation over the block + endpoints
            block = linspace(v_before, v_after, len+2);
            % block is 1×(len+2), indexes:
            % block(1) is v_before, block(end) is v_after
            % the middle len entries correspond to the NaNs
            y(i0:i1) = block(2:end-1);
        end
    end

    if wasRow
        y = y.';
    end
end
