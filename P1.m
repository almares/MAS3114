diary on
format compact

%Exercise 1:
%1:
randi(2,5,4)-1
ans =
     0     1     0     0
     1     0     1     1
     0     0     0     1
     1     0     0     1
     1     1     0     0
randi(2,5,4)-1
ans =
     1     0     0     0
     1     0     0     0
     0     0     1     1
     1     0     1     0
     0     0     0     0
randi(2,5,4)-1
ans =
     1     0     0     0
     0     1     0     1
     0     1     1     1
     0     1     0     1
     0     0     0     0

%2
x=randi(100,4,1)
x =
    18
    23
    44
    32
X=[x.^0 x.^1 x.^2 x.^3 x.^4 x.^5]
X =
  Columns 1 through 4
           1          18         324        5832
           1          23         529       12167
           1          44        1936       85184
           1          32        1024       32768
  Columns 5 through 6
      104976     1889568
      279841     6436343
     3748096   164916224
     1048576    33554432

%3
j=flipud(diag((randi(10,5,1)-1)))
j =
     0     0     0     0     4
     0     0     0     1     0
     0     0     8     0     0
     0     4     0     0     0
     9     0     0     0     0

%4
J=diag(randi(90,6,1)+10)
J =
    86     0     0     0     0     0
     0    40     0     0     0     0
     0     0    60     0     0     0
     0     0     0    99     0     0
     0     0     0     0    60     0
     0     0     0     0     0    40
J(1,:)=randi(90,1,6)+10
J =
    66    43    79    48    55    73
     0    40     0     0     0     0
     0     0    60     0     0     0
     0     0     0    99     0     0
     0     0     0     0    60     0
     0     0     0     0     0    40
J(:,1)=randi(90,6,1)+10
J =
    98    43    79    48    55    73
    40    40     0     0     0     0
    86     0    60     0     0     0
    77     0     0    99     0     0
    96     0     0     0    60     0
    13     0     0     0     0    40

%Exercise 2:
type multi

function [C, CRows, CColumns] = multi(A,B)
%This function finds the product of 2 matricies using
%the definition of product, row-by-row, and column-ny
%column methods for calculation. If the dimensions are
%inconsistent, it returns an error

[m,p] = size(A);
[q,n] = size(B);

C=zeros(m,n);
CRows=zeros (m,n);
CColumns=zeros(m,n);

if p==q,
    %definition
    for i=1 : n
        a=B(:,i);
        b=A*a;
        C(:,i)=b;
    end
    %row-by-row
    for i=1 : m
        for j=1 : n
            temp = 0;
            for k=1 : p
                temp = temp + A(i,k)*B(k,j);
                CRows(i,j) = temp; 
            end
        end
    end  
    %column-by-column
    for i=1 : n
        for j=1 : m
            temp = 0;
            for k=1 : p
                temp = temp + A(j,k)*B(k,i);
                CColumns(j,i) = temp; 
            end
        end
    end  
else
    C=[];
    CRows=[];
    CColumns=[];
    disp('The dimensions of A and B disagree')
end
end

%(a)
A=randi(10,2,3)
A =
     9     2     7
    10    10     1
B=magic(2)
B =
     1     3
     4     2
[C,CRows,CColumns]= multi(A,B)
The dimensions of A and B disagree
C =
     []
CRows =
     []
CColumns =
     []
%(b)
A= magic(5) 
A =
    17    24     1     8    15
    23     5     7    14    16
     4     6    13    20    22
    10    12    19    21     3
    11    18    25     2     9
B= ones(4,6)
B =
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
[C,CRows,CColumns]= multi(A,B)
The dimensions of A and B disagree
C =
     []
CRows =
     []
CColumns =
     []
&(c)
A = magic(4)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
B = ones(4,3)
B =
     1     1     1
     1     1     1
     1     1     1
     1     1     1
[C,CRows,CColumns]= multi(A,B)
C =
    34    34    34
    34    34    34
    34    34    34
    34    34    34
CRows =
    34    34    34
    34    34    34
    34    34    34
    34    34    34
CColumns =
    34    34    34
    34    34    34
    34    34    34
    34    34    34
%(d)
A = ones(4)
A =
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
B = diag([2,3,4,5])
B =
     2     0     0     0
     0     3     0     0
     0     0     4     0
     0     0     0     5
[C,CRows,CColumns]= multi(A,B)
C =
     2     3     4     5
     2     3     4     5
     2     3     4     5
     2     3     4     5
CRows =
     2     3     4     5
     2     3     4     5
     2     3     4     5
     2     3     4     5
CColumns =
     2     3     4     5
     2     3     4     5
     2     3     4     5
     2     3     4     5

A=randi(10,2,3),B=magic(2) 
A =
     3    10     2
     6    10    10
B =
     1     3
     4     2
A*B
{_Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('mtimes')" style="font-weight:bold"> * </a>
Inner matrix dimensions must agree.}_ 

A= magic(5), B= ones(4,6)
A =
    17    24     1     8    15
    23     5     7    14    16
     4     6    13    20    22
    10    12    19    21     3
    11    18    25     2     9
B =
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
     1     1     1     1     1     1
A*B
{_Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('mtimes')" style="font-weight:bold"> * </a>
Inner matrix dimensions must agree.}_ 

A = magic(4), B = ones(4,3)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
B =
     1     1     1
     1     1     1
     1     1     1
     1     1     1
A*B
ans =
    34    34    34
    34    34    34
    34    34    34
    34    34    34

A = ones(4), B = diag([2,3,4,5])
A =
     1     1     1     1
     1     1     1     1
     1     1     1     1
     1     1     1     1
B =
     2     0     0     0
     0     3     0     0
     0     0     4     0
     0     0     0     5
A*B
ans =
     2     3     4     5
     2     3     4     5
     2     3     4     5
     2     3     4     5

%The outputs from multi(A,B) match the multiplication function in MATLAB. The errors from both functions are identical.
diary off

diary on
%Exercise 3
type givensrot

function G = givensrot(n,i,j,theta)
%This function produces an n*m Givens matrix. Returns an empty matrix if
%conditions not met and a Givens rotation matrix otherwise.

G = eye(n);

if 1 <= i && i < j && i <= n
    G(i,i) = cos(theta);
    G(j,j) = cos(theta);
    G(j,i) = sin(theta);
    G(i,j) = -sin(theta);
    return,
end

if j>i
    G = [];
    return,
end
end

%a1)
givensrot(4,3,2,pi/2)

ans =

     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
 
%a2)
givensrot(5,2,4,pi/4)

ans =

    1.0000         0         0         0         0
         0    0.7071         0   -0.7071         0
         0         0    1.0000         0         0
         0    0.7071         0    0.7071         0
         0         0         0         0    1.0000
 
%a3)
givensrot(3,1,2,pi)

ans =

   -1.0000   -0.0000         0
    0.0000   -1.0000         0
         0         0    1.0000
 
%b)
e1=[1;0;0]

e1 =

     1
     0
     0

givensrot(3,1,2,pi)*e1

ans =

   -1.0000
    0.0000
         0

e2=[0;1;0]

e2 =

     0
     1
     0

givensrot(3,1,2,pi)*e2

ans =

   -0.0000
   -1.0000
         0

e3=[0;0;1]

e3 =

     0
     0
     1

givensrot(3,1,2,pi)*e3

ans =

     0
     0
     1
 
x=[1;1;1]

x =

     1
     1
     1

givensrot(3,1,2,pi)*x

ans =

   -1.0000
   -1.0000
    1.0000
 
%These singular column outputs are consistent with the geometrical meaning, indicating a rotation in the x-plane exclusively through theta radians. When multiplied by each vector, the matrix G was only affected in terms of the x variable in the initial column.

diary off

% Exercise 4
diary on
format compact
type toeplitz
function A=toeplitz(m,n,a)
A=zeros(m,n);
    for i=1:m
        for j=1:n
            A(i,j)=a(nj);
        end
    end
end
%(1)
%(a) 
m=4;n=3;a=transpose([1:6])
a =
     1
     2
     3
     4
     5
     6
A=toeplitz(m,n,a)
A =
     3     2     1
     4     3     2
     5     4     3
     6     5     4
%(b) 
m=3;n=4;a=randi(10,6,1)
a =
     3
     6
    10
    10
     2
    10
A=toeplitz(m,n,a)
A =
    10    10     6     3
     2    10    10     6
    10     2    10    10
%(c)
m=4;n=4;a=[zeros(3,1);[1:4]']
a =
     0
     0
     0
     1
     2
     3
     4
A=toeplitz(m,n,a)
A =
     1     0     0     0
     2     1     0     0
     3     2     1     0
     4     3     2     1
%(2)
m=5;n=5;a=[randi([10,100],5,1);zeros(4,1)]
a =
    84
    32
    94
    41
    27
     0
     0
     0
     0
A=toeplitz(m,n,a)
A =
    27    41    94    32    84
     0    27    41    94    32
     0     0    27    41    94
     0     0     0    27    41
     0     0     0     0    27

%(3)
a=[0 0 0 1 0 0 0]'
a =
     0
     0
     0
     1
     0
     0
     0
m=4;n=4
n =
     4
A=toeplitz(m,n,a)
A =
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
diary off

diary on
format compact

%Exercise#5

type stochastic
 
function P=stochastic(A)
    a = false;
    b = false;
    [m,n] = size(A);
    for i=1:m
        y = zeros(m);
        z = zeros(n);
        for j=1:n
            y(i,j) = A(i,j);
            z(i,j) = A(j,i);
        end
        if y==zeros(m)
            a = true;
        end
        if z == zeros(n)
            b = true;
        end
	  %loops check if the square matrix contains a zero column and zero row
    end
    if(a==true && b== true)
        disp('A is not stochastic and cannot be scaled to stochastic')
        P=[];
    end
    c = false;
    d = false;
    e = false;
    %Starting with the assumption that the rows and columns are probability             vectors
    for i=1:m
        total = sum(A(i,:));
        if total~=1
            c=true;
        end
    end
    for j=1:m
        total = sum(A(:,j));
        if total~=1
            d=true;
        end
    end
    if c==true
        if d==false
            disp('It is left stochastic’);
		% Only columns are probability vectors
            P=A;
        elseif d==true && (a==false || b==false)
		disp('A is not stochastic and can be scaled to stochastic')
            if b==true
            for i=1:m
                total = sum(A(i,:));
                multip = (1/total);
                for j=1:n
                    P(i,j)=(multip*(A(i,j)));
                end
            end
            elseif a==true
                for j=1:n
                total = sum(A(:,j));
                multip = (1/total);
                for i=1:n
                    P(i,j)=(multip*(A(i,j)));
                end
                end
            elseif a==false && b==false
                for j=1:n
                total = sum(A(:,j));
                multip = (1/total);
                for i=1:n
                    P(i,j)=(multip*(A(i,j)));
                end
                end
            end
        end
    elseif d==true
        disp('It is right stochastic’);
	  % Only rows are probability vectors
            P=A;
    elseif c==false && d==false
        disp('It is doubly stochastic');
            P=A;
	  % Both Rows and Columns are probability vectors
    end
end

% part a
A=[0.5,  0,  0.5; 0,  0,  1; 0.5,  0,  0.5]  
A =
    0.5000         0    0.5000
         0         0    1.0000
    0.5000         0    0.5000
P=stochastic(A)
It is right stochastic
P =
    0.5000         0    0.5000
         0         0    1.0000
    0.5000         0    0.5000
% part b
A=A'
A =
    0.5000         0    0.5000
         0         0         0
    0.5000    1.0000    0.5000
P=stochastic(A)
It is left stochastic
P =
    0.5000         0    0.5000
         0         0         0
    0.5000    1.0000    0.5000

% part c
A=[0.5,  0,  0.5; 0,  0,  1; 0,  0,  0.5]  
A =
    0.5000         0    0.5000
         0         0    1.0000
         0         0    0.5000
P=stochastic(A)
P =
    0.5000         0    0.5000
         0         0    1.0000
         0         0    1.0000
% part d
A=A'
A =
    0.5000         0         0
         0         0         0
    0.5000    1.0000    0.5000
P=stochastic(A)
P =
    0.5000         0         0
         0         0         0
    0.5000    1.0000    1.0000
% part e
A=[0.5,  0,  0.5; 0,  0.5,  0.5; 0.5,  0.5,  0]  
A =
    0.5000         0    0.5000
         0    0.5000    0.5000
    0.5000    0.5000         0
P=stochastic(A)
It is doubly stochastic
P =
    0.5000         0    0.5000
         0    0.5000    0.5000
    0.5000    0.5000         0

% part f
A=magic(3)
A =
     8     1     6
     3     5     7
     4     9     2
P=stochastic(A)
P =
    0.5333    0.0667    0.4000
    0.2000    0.3333    0.4667
    0.2667    0.6000    0.1333
% part g
A= diag([1,2,3])
A =
     1     0     0
     0     2     0
     0     0     3
P=stochastic(A)
P =
     1     0     0
     0     1     0
     0     0     1

%part g
A=[0,  0,  0; 0,  0.5,  0.5; 0,  0.5,  0.5]  
A =
         0         0         0
         0    0.5000    0.5000
         0    0.5000    0.5000
P=stochastic(A)
A is not stochastic and cannot be scaled to stochastic
P =[]

diary off

diary Project1
diary on
format compact

%Exercise6

x=linspace(0,4,8);
y=atan(x)+x-1;
 
plot(x,y) 
x=0.5;
syms x
f=atan(x)+x-1

f =
x + atan(x) – 1

g=diff(f)

g =
1/(x^2 + 1) + 1

type newtons
function root=newtons(N,x)
    for n=1:N
        x1 = (x-(atan(x)+x-1)/(1/(x^2+1)+1));
        x=x1;
        fprintf('%.12f\n',x1)
    end
    root = x1;
end

N=5;
x=0.5;

root=newtons(N,x)

0.520195772777
0.520268991753
0.520268992720
0.520268992720
0.520268992720
root =
    0.5203

% The intermediate outputs were all very close to the actual root and increased steeply in accuracy
%The approximation of this function is 0.52026899
x=linspace(0,4,8);
y=x.^3-x-1;
 
plot(x,y) 

syms x
f=x^3-x-1

f =
x^3 - x – 1

g=diff(f)

g =
3*x^2 – 1

format compact
N=5;
x=1.5;

root=newtons_1(N,x)

1.347826086957
1.325200398951
1.324718173999
1.324717957245
1.324717957245
root =
    1.3247

%The outputs took a little longer to zero in on a solution
%The approximation seems to be 1.32471795
x=1;

root=newtons_1(N,x)

1.500000000000
1.347826086957
1.325200398951
1.324718173999
1.324717957245
root =
    1.3247

%It's interesting that the first value was exactly 1.5 and then continued as had been done in the previous test
%The approximation seems to still be 1.32471795
x=0.6;

root=newtons_1(N,x)

17.900000000000
11.946802328609
7.985520351936
5.356909314795
3.624996032946
root =
    3.6250

%The values are very skewed now and the variations are very large
%The approximation based on this data is incorrect but it would be 3.62499603
x=0.57;

root=newtons_1(N,x)

-54.165454545454
-36.114292524926
-24.082094252098
-16.063387407818
-10.721483416797
root =
  -10.7215

%The values are way off now
%Based off this skewed information the approximation is -10.72148341
N=10;
x=0.6;

root=newtons_1(N,x)

17.900000000000
11.946802328609
7.985520351936
5.356909314795
3.624996032946
2.505589190107
1.820129422319
1.461044109888
1.339323224263
1.324912867719
root =
    1.3249

x=0.57;

root=newtons_1(N,x)

-54.165454545454
-36.114292524926
-24.082094252098
-16.063387407818
-10.721483416797
-7.165534466882
-4.801703812713
-3.233425234527
-2.193674204845
-1.496866569238
root =
   -1.4969

N=100;
x=0.6;

root=newtons_1(N,x)

1.324717957245
root =
    1.3247

x=0.57;

root=newtons_1(N,x)

1.324717957245
root =
    1.3247

%The initial process does eventually converge yet it is a very slow process for part 4 as opposed to part 3
% The reason is that less than 0.577, the graph is decreasing and therefore the derivative is negative, messing up the equation.
%So anything above 0.577 is fast yet everything below it takes a long time. 

diary off


