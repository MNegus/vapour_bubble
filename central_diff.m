function diff = central_diff(arr, dx)
%CENTRAL_DIFF Returns the central difference array
%   Function to return the central difference array, an approximation to
%   the first derivative. If arr is length N, then the array diff is of
%   length N - 2, as the central difference isn't defined at the end
%   points.
    
N = length(arr);
diff = (arr(3 : N) - arr(1 : N - 2)) / (2 * dx);
end

