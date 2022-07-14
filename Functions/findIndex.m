function closestIndex = findIndex(matrix, numberToFind)
if numberToFind > 1
    closestIndex = zeros(1, numel(numberToFind));
    for i = 1 : length(numberToFind)
        [~, closestIndex(i)] = min(abs(matrix - numberToFind(i)));
    end
else
    [~, closestIndex] = min(abs(matrix - numberToFind));
end
end