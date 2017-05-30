function y = nanstd( x, dim, flag )

% Compute the standard deviation while ignoring NaN values
% If all the values are NaN, the standard deviation is returned as NaN
% If there is only a single non-NaN value, the deviation is returned as 0

% INPUTS:
%   x    : array with inout data
%   dim  : Obtaining the standard deviation along the dimension dim of x
%   flag : If flag=1 output values will be normilized by n
%          If flag=0 output values are normilizaed by n-1
%          (n=number of remaining observations)

% Check inputs 
if isempty(x)
	y = NaN;
	return
end

if nargin < 3
	flag = 0;
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1; 
	end	  
end

% Find NaNs in x and nanmean(x)
nans = isnan(x);
avg = nanmean(x,dim);

% create array indicating number of element 
% of x in dimension DIM (needed for subtraction of mean)
tile = ones(1,max(ndims(x),dim));
tile(dim) = size(x,dim);

% remove mean
x = x - repmat(avg,tile);

count = size(x,dim) - sum(nans,dim);

% Replace NaNs with zeros.
x(isnan(x)) = 0; 

% Protect against a  all NaNs in one dimension
i = find(count==0);

if flag == 0
	y = sqrt(sum(x.*x,dim)./max(count-1,1));
else
	y = sqrt(sum(x.*x,dim)./max(count,1));
end
y(i) = i + NaN;