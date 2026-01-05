function A = k_regular_graph(N, k, maxSwaps)
if mod(k*N,2) ~= 0
    error('k*N must be even.');
end
if k >= N
    error('k must be less than N.');
end
if mod(k,2) ~= 0
    error('This version supports only even k.');
end
% Step 1: Generate initial k-regular graph
A = zeros(N);
for i = 1:N
    for j = 1:(k/2)
        neighbor = mod(i + j - 1, N) + 1;
        A(i, neighbor) = 1;
        A(neighbor, i) = 1;
    end
end
% Step 2: Randomize with edge swaps
edgeList = find(triu(A)); % upper triangle indices
[r, c] = ind2sub([N N], edgeList);
E = [r c];
swaps = 0;
attempts = 0;
while swaps < maxSwaps && attempts < 10 * maxSwaps
    attempts = attempts + 1;
    % Pick two random edges
    idx = randperm(size(E, 1), 2);
    [a, b] = deal(E(idx(1), 1), E(idx(1), 2));
    [c, d] = deal(E(idx(2), 1), E(idx(2), 2));
    % Make sure all nodes are distinct
    if length(unique([a, b, c, d])) < 4
        continue;
    end
    % Proposed new edges
    if rand < 0.5
        u1 = a; v1 = d;
        u2 = c; v2 = b;
    else
        u1 = a; v1 = c;
        u2 = d; v2 = b;
    end
    % Skip if any proposed edge already exists or is a self-loop
    if A(u1, v1) || A(u2, v2) || u1 == v1 || u2 == v2
        continue;
    end
    % Perform the swap
    A(a, b) = 0; A(b, a) = 0;
    A(c, d) = 0; A(d, c) = 0;
    A(u1, v1) = 1; A(v1, u1) = 1;
    A(u2, v2) = 1; A(v2, u2) = 1;
    E(idx(1), :) = sort([u1, v1]);
    E(idx(2), :) = sort([u2, v2]);
    swaps = swaps + 1;
end
