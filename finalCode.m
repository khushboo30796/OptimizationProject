original = double(imread('example.bmp'))/255;
marked   = double(imread('example_marked.bmp'))/255;

annotate_thresh = 0.25; % may need to be adjusted
colorim = sum(abs(original - marked), 3) > annotate_thresh;
colorim = double(colorim);
spy(colorim);
anchorind = find(colorim);

[m,n,p] = size(marked);
assert(p == 3);
N = m*n;

% convert from RGB to YUV using the matrix T
T = [1.0 0.956 0.621; 1.0 -0.272 -0.647; 1.0 -1.106 1.703];
yuv = reshape(marked(:), N, 3)/T';
yuv = reshape(yuv, m, n, 3);
chromi = yuv(:,:,2);
chromq = yuv(:,:,3);
% Get grayscale value from original image
yuv = reshape(original(:), N, 3)/T';
yuv = reshape(yuv, m, n, 3);
gray = yuv(:,:,1);

% Get anchors from the I and Q color components
anchor(:,1) = chromi(anchorind);
anchor(:,2) = chromq(anchorind);

edges=[(1:N)' ((1:N)+1)'];        %add downward edge between node i with node i+1
edges=[edges; (1:N)' (1:N)'+m];   %add leftward edge between node i and i+m
% exclude edges which link to nodes outside the domain
excluded=edges(:,1)>N | edges(:,1)<1 | edges(:,2)>N |edges(:,2)<1; 
edges(excluded, :) = [];
% remove edge from bottom of column to top of next column
edges((m:m:(m-1)*n)',:)=[]; 

beta = 200;
EPSILON = 1e-5;

% calculate edge weights
graydistance = abs(gray(edges(:,1)) - gray(edges(:,2)));
ming = min(graydistance);  % it's helpful to normalize distances
maxg = max(graydistance);
graydistance = (graydistance - ming)./(maxg - ming);
weights=exp(-(beta*graydistance)) + EPSILON;

%Build sparse weighted adjacency matrix
W=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)], ...
    [weights;weights],N,N);
% Build Degree matrix
D = diag(sum(W));
% Build Laplacian
L=D-W;

mark = anchorind;
rest = 1:n*m;
rest(mark) = [];
rest = rest';
% Extract the portion of the Lapalcian matrices
Lu = L(rest, rest);
R  = L(rest, mark);
recover = zeros(m, n, 3);
recover(:,:,1) = gray;
xchannel = zeros(m*n,1);
for i=1:2
    d  = anchor(:,i);
    xrest = -Lu\(R*d);
    xchannel(rest) = xrest;
    xchannel(mark) = d;
    recover(:,:,i+1) = reshape(xchannel, m, n);
end
% Convert from YUV back to RGB using T matrix
recover = reshape(recover, m*n, 3);
recover = recover*T';
% Keep RGB values within [0,1]
recover = min(max(0, recover), 1);
final   = reshape(recover, m, n, 3);

figure(1);
image(final);