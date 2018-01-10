function x = vectorize(X)
% Reshapes a matrix into a vector
x = reshape(X, [numel(X), 1]);
end