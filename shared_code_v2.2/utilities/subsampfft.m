function Y = subsampfft(X,S,mode,im_sz)


% S is the sampled coefficients
% mode = 1 for forward operator
% mode = 2 for conjugate transpose
% im_sz is the original image size

switch mode
    case 1
        X = reshape(X,im_sz);
        Y = fftn(X);%/sqrt(numel(X));
        Y = Y(S);
        Y = Y(:);
        
    case 2
        Z = zeros(im_sz);
        Z(S) = X;
        Y = ifftn(Z)*numel(Z);%sqrt(numel(Z));
        Y = Y(:);
        
end