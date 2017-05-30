function means = nanmean( data, dim )

% This function calculates the mean of a data containing NaN values
%
% INPUTS:
%   data : It is the array containing the data
%   dim  : 1=mean value for each column; 2=mean value for each row
%
% OUTPUT:
%   means : mean value

% If there is only 1 input, it is assumed array [nx1]
if nargin == 1
    dim = 1;
end

if dim == 1 %(columns)
    dim_length = length(data(1,:));
    means = zeros(1,dim_length);
    for i = 1:dim_length
        col = data(:,i);
        m = mean(col(~isnan(col)));
        means(i) = m;
    end
elseif dim == 2 %(rows)
    dim_length = length(data(:,1));
    means = zeros(dim_length,1);
    for i = 1:dim_length
        col = data(i,:);
        m = mean(col(~isnan(col)));
        means(i) = m;
    end
else
    disp('Error calculating nanmean')
end
