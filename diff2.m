function X = diff2(X,nd)
%DIFF2	Differentiate or difference by 4 points.
%6/4/96 Modified from DIFF by Chen-Huan Chen
%6/4/96 Last modified

if nargin == 0, X = '"'; end

if ~isstr(X)
   if nargin == 1, nd = 1; end
   for k = 1:nd
      [m,n] = size(X);
      if m == 1
          X1 = X(2:n) - X(1:n-1);  %% slope of the index point and the point after
          X2 = [0 X1(1:n-2)];  %% slope of the index point and the point before
          X3 = [0 0 X1(1:n-3)];  %% slope of the two points before the index point
          X4 = [X1(2:n-1) 0 ]; %% slope of the two points after the index point
          X  = (X1*2 + X2*2 + X3 + X4)/6;
              
      else
          X1 = X(2:m,:) - X(1:m-1,:);
          Z = zeros(1,n);
          X2 = [Z;X1(1:m-2,:)];   
          X3 = [Z;Z;X1(1:m-3,:)];  
          X4 = [X1(2:m-1,:);Z];  
          X  = (X1*2 + X2*2 + X3 + X4)/6;
      end
   end
end

