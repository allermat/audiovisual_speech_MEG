function [cxy,fc] = coherence(x,y,fs,varargin)
% Compute magnitude-squared coherence with optional partialization
% 
% Input: 
%   x, y: The two (sets of) signals between which to compute coherence. 
%         Number of samples along first dimension, number of signals along
%         second dimension. Can be a single column vector. If both are
%         matrices, sizes must match. 
%   fs: sampling frequency (scalar)
%   partialize: The signal to be partialized out. 
%   square: Wheter to compute the coherence squaring the magnitude or not.
%       Deafault is square = true
% Output: 
%   cxy: the magnitude-square coherence spectrum between x and y
%   fc: vector of frequencies at which coherence is evaluated
% 
% Copyright(C) Mate Aller 2020
% allermat@gmail.com

p = inputParser;

addRequired(p,'x',@(x) validateattributes(x,{'numeric'},{'2d'}));
addRequired(p,'y',@(x) validateattributes(x,{'numeric'},{'2d'}));
addRequired(p,'fs',@(x) validateattributes(x,{'numeric'},{'scalar','positive'}));
addParameter(p,'partialize',[],@(x) validateattributes(x,{'numeric'},{'2d'}));
addParameter(p,'square',true,@(x) validateattributes(x,{'logical'},{'scalar'}));
parse(p,x,y,fs,varargin{:});

x = p.Results.x;
y = p.Results.y;
fs = p.Results.fs;
z = p.Results.partialize;
square = p.Results.square;

% The number of samples in all input signals (first dimension) must be equal
nSample = size(x,1);
assert(size(y,1) == nSample)
if ~isempty(z)
    assert(size(z,1) == nSample);
end

[pxy,fc] = cpsd(x,y,[],[],nSample,fs);
pxx = cpsd(x,x,[],[],nSample,fs);
pyy = cpsd(y,y,[],[],nSample,fs);

% Optionally partialize the signal z out
if ~isempty(z)
    pxz = cpsd(x,z,[],[],nSample,fs);
    pyz = cpsd(y,z,[],[],nSample,fs);
    pzz = cpsd(z,z,[],[],nSample,fs);

    [pxx,pxy,pyy] = partialize_csd(pxx,pxy,pxz,pyy,pyz,pzz);
end
% Compute magnitude-squared coherence
if square
    cxy = abs(mean(pxy,2)).^2./(abs(mean(pxx,2)).*abs(mean(pyy,2)));
else
    cxy = abs(mean(pxy,2))./sqrt((abs(mean(pxx,2)).*abs(mean(pyy,2))));
end
end

function [pxx_p,pxy_p,pyy_p] = partialize_csd(pxx,pxy,pxz,pyy,pyz,pzz)
% partial spectra are computed as in Rosenberg JR et al (1998) 
% J.Neuroscience Methods, equation 38

% size of cross-power spectra: nfreq x nrpt
siz = size(pxy);

% The signal in z will be partialled out
A = zeros(siz(2),3,3,siz(1));
for j = 1:siz(2) % loop over rpt
    AA = cat(1,cat(2,shiftdim(pxx(:,j),-2),shiftdim(pxy(:,j),-2),shiftdim(pxz(:,j),-2)),...
               cat(2,shiftdim(pxy(:,j),-2),shiftdim(pyy(:,j),-2),shiftdim(pyz(:,j),-2)),...
               cat(2,shiftdim(pxz(:,j),-2),shiftdim(pyz(:,j),-2),shiftdim(pzz(:,j),-2)));
    AB = cat(1,shiftdim(pxz(:,j),-2),shiftdim(pyz(:,j),-2),shiftdim(pzz(:,j),-2));
    BA = cat(2,shiftdim(pxz(:,j),-2),shiftdim(pyz(:,j),-2),shiftdim(pzz(:,j),-2));
    BB = shiftdim(pzz(:,j),-2);
    for k = 1:siz(1) % loop over freq
        A(j,:,:,k) = AA(:,:,k) - AB(:,:,k)/(BB(:,:,k))*BA(:,:,k);
    end
end

pxx_p = squeeze(A(:,1,1,:))';
pxy_p = squeeze(A(:,2,1,:))';
pyy_p = squeeze(A(:,2,2,:))';

end