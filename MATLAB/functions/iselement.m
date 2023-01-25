function result = iselement(M, condition)

    valid = zeros(size(M));
    for c = 1:length(condition)
        valid = valid | (M == condition(c));
    end

    result = valid;
end