% X: denoised matrix, s2: original noise variance, p: number of signal components, s2_after: noise variance after denoising
function [X,s2,p,s2_after] = denoiseMatrix(X) 
M = size(X,1);
N = size(X,2);
if M<N
    [U,lambda] = eig(X*X','vector');
else
    [U,lambda] = eig(X'*X,'vector');
end
[lambda,order] = sort(lambda,'descend');
U = U(:,order);
csum = cumsum(lambda,'reverse');
p = (0:length(lambda)-1)';
p = -1 + find((lambda-lambda(end)).*(M-p).*(N-p) < 4*csum*sqrt(M*N),1);
if p==0
    X = zeros(size(X));
elseif M<N
    X = U(:,1:p)*U(:,1:p)'*X;
else
    X = X*U(:,1:p)*U(:,1:p)';
end
s2 = csum(p+1)/((M-p)*(N-p));
s2_after = s2 - csum(p+1)/(M*N);
end