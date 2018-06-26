function diff = second_diff(arr,dx)
%SECOND_DIFF Returns the second difference array
%   Function to return the second difference array, an approximation to
%   the second derivative. If arr is length N, then the array diff is of
%   length N - 2, as the second difference isn't defined at the end points.

N = length(arr);
diff = (arr(3:N) - 2 * arr(2:N-1) + arr(1:N-2)) / (dx^2);
end

