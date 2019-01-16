function	s=cumsimpsum(y,DEL)
% CUMSIMPSUM Cumulative numerical integration via the Simpson's rule.
% 	Z = CUMSIMPSUM(Y,DEL) computes an approximation of the cumulative
% 	integral of Y from a=0 to b via the Simpson's rule (with unit spacing).
% 	For matrices Y, CUMSIMPSUM works on columns of Y. CUMSIMPSUM needs at
% 	least 3 elements in the column(s) of Y. DEL contains the upper limit(s)
% 	b. The default value of the upper limit equals 1. DEL is a scalar if the
% 	upper limit for all column vectors in Y is the same. DEL is a row vector
% 	if the upper limits for the column vectors in Y are different. CUMSIMPSUM
% 	assumes the remaining upper limits to equal 1 if the elements in DEL are
% 	not equal to the number of columns of Y.
% 	
% 	The area of the first panel is calculated based on an parabolic
% 	estimation. The area of the remaining panels is calculated based on
% 	Simpson's 1/3 and 3/8 Rule. The cumulative sum of the panels is
% 	calculated.
% 	
% 	Example: If 
% 	
% 	N=8;
% 	X=ones(1,N);
% 	Y=[X',X'];
% 	DEL=[N-1,2*(N-1)];
% 	
% 	then
% 	
% 	S=cumsimpsum(Y,DEL);
%                
% 	S =     0     0
%           1     2
%           2     4
%           3     6
%           4     8
%           5    10
%           6    12
%           7    14
% 	
% 	See also SIMPSUM in MATLAB FILE EXCHANGE, CUMTRAPZ.
% 	
% 	Michalis Meyer, Oct. 28th 2004
% 	Thanks to Vassili Pastushenko for his suggestions.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                       START OF EXECUTABLE CODE                      %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b]=size(y);
if a==1
    y=y(:);
end

[N,b]=size(y);

if N<3
    warning('CUMSIMPSUM needs at least 3 elements in the column of Y.');
end

% Upper limits
if nargin==1
    DEL=ones(1,b);
elseif nargin==2
    [c,d]=size(DEL);
	if (c==1)
		if d==1
            DEL=DEL*ones(1,b);
		elseif (d~=1)&(d<b)
            DEL=[DEL,ones(1,b-d)];
		elseif (d~=1)&(d>b)
            warning('Number of elements in row vector DEL is larger than number of columns in Y');
		end;
	else warning('DEL needs to be a scalar or row vector.');
	end
end

% Spacing
h=DEL/(N-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                                                     %
%                         Numerical Integration                       %
%                                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% First two panels
s=zeros(size(y));
s(2,:)=h/12.*(5*y(1,:)+8*y(2,:)-y(3,:));
s(3,:)=h/3.*(y(1,:)+4*y(2,:)+y(3,:));

% Remaining panels
if N>3
	% Computing approximation of the cumulative integral
	if rem(N,2)==0 % N even, odd # of panels
        
        s1=(y(3:2:N-3,:)+4*y(4:2:N-2,:)+y(5:2:N-1,:));s1=ones(size(s1,1),1)*h/3.*cumsum(s1,1);
        s2=(y(4:2:N-2,:)+4*y(5:2:N-1,:)+y(6:2:N,:));s2=ones(size(s2,1),1)*h/3.*cumsum(s2,1);
        
        s(4,:)=3*h/8.*(y(1,:)+3*y(2,:)+3*y(3,:)+y(4,:));
        s(5:2:N-1,:)=s1+ones(size(s1,1),1)*s(3,:);
        s(6:2:N,:)=s2+ones(size(s2,1),1)*s(4,:);
        
	else % N odd, even # of panels
        
        s1=(y(3:2:N-2,:)+4*y(4:2:N-1,:)+y(5:2:N,:));s1=ones(size(s1,1),1)*h/3.*cumsum(s1,1);
        s2=(y(4:2:N-3,:)+4*y(5:2:N-2,:)+y(6:2:N-1,:));s2=ones(size(s2,1),1)*h/3.*cumsum(s2,1);
        
        s(4,:)=3*h/8.*(y(1,:)+3*y(2,:)+3*y(3,:)+y(4,:));
        s(5:2:N,:)=s1+ones(size(s1,1),1)*s(3,:);
        s(6:2:N-1,:)=s2+ones(size(s2,1),1)*s(4,:);
        
	end
end

if a==1
    s=s';
end