function laplacian = lap3d(n)
% Generate discrete 3D laplacian matrix for a flattened real-space vector
% with n elements. Returns an n^3xn^3 matrix.

blockA = diag(ones(1,n)*-6.) + diag(ones(1,n-1),1) + diag(ones(1,n-1),-1);
Ar = repmat(blockA, 1, n); % repeat matrix for block diagonal
Ac = mat2cell(Ar, size(blockA,1), repmat(size(blockA,2),1,n));
blockB = blkdiag(Ac{:}) + diag(ones(1,n^2-n),n) + diag(ones(1,n^2-n),-1*n);

Br = repmat(blockB,1,n); % repeat matrix for block diagonal
Bc = mat2cell(Br, size(blockB,1), repmat(size(blockB,2),1,n));
laplacian = blkdiag(Bc{:}) + diag(ones(1,n^3-n^2),n^2) +diag(ones(1,n^3-n^2),-1*n^2);
