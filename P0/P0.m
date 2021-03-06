diary on
format compact

%Exercise 1
diary off


%Exercise 2
A = [-2 -1 8 3
0 -2 3 1
-3 -7 5 4] 

A =
    -2    -1     8     3
     0    -2     3     1
    -3    -7     5     4

B = [1 -2
3 -1
1 0
]

B =
     1    -2
     3    -1
     1     0

X = [-2
3
5]

X =
    -2
     3
     5

x = [1
2
3
4]

x =
     1
     2
     3
     4

x = [1 2 3 4]

x =
     1     2     3     4

y = [1
4
6
3]

y =
     1
     4
     6
     3

size(A)
ans =
     3     4
size(B)
ans =
     3     2
size(X)
ans =
     3     1
size(x)
ans =
     1     4
size(y)
ans =
     4     1
%The matrices x and y are not of the same size, as their number of rows and columns differ. The command size() of x and y returned the same numbers, but in different orders.

F = [5 2 -3
4 3 -2
0 -1 6
1 0 -2]
F =
     5     2    -3
     4     3    -2
     0    -1     6
     1     0    -2

A(1,3)
ans =
     8
%The command A(1,3) returned the element in row one and column 3: 8
A(:,3)
ans =
     8
     3
     5
%The command A(:,3) returned the elements in column 3 of matrix A 
A(2,:)
ans =
     0    -2     3     1
%The command A(2,:) returned the elements in row 2 of matrix A


A([1 2],[3 4])
ans =
     8     3
     3     1
%The command A([1 2], [3 4]) returned the top left corner four elements of matrix A
F(:,4)=[-1 1 -4 3]
F =
     5     2    -3    -1
     4     3    -2     1
     0    -1     6    -4
     1     0    -2     3
%The command F(:,4)=[-1 1 -4 3] returned an altered matrix of F with an added column with the four elements entered in
F([1 3],[2 4])=[1 -3; 2 -4]
F =
     5     1    -3    -3
     4     3    -2     1
     0     2     6    -4
     1     0    -2     3
%The command F([1 3],[2 4])=[1 -3; 2 -4] replaced the elements in rows one and three and columns two and four with 1, -3, 2, and -4 respectively.
F([2 3], :)=A([1 3],:)
F =
     5     1    -3    -3
    -2    -1     8     3
    -3    -7     5     4
     1     0    -2     3
%The command F([2 3], :)=A([1 3],:) returned an altered F matrix, replacing its second and third rows with the second and third rows of matrix A.
F(:,[1 2])=F(:,[2 1])
F =
     1     5    -3    -3
    -1    -2     8     3
    -7    -3     5     4
     0     1    -2     3
%The command F(:,[1 2])=F(:,(2 1]) returned an altered F matrix, switching its first and second columns.

F(:,[1])=y
F =
     1     5    -3    -3
     4    -2     8     3
     6    -3     5     4
     3     1    -2     3
%Replaces first column of matrix F with vector y.

F([1 2],:)=F([2 1],:)
F =
     4    -2     8     3
     1     5    -3    -3
     6    -3     5     4
     3     1    -2     3
%Switches rows one and two of matrix F.

[A B]
ans =
    -2    -1     8     3     1    -2
     0    -2     3     1     3    -1
    -3    -7     5     4     1     0
%This combined matrices A and B together, attaching B to the end of A.

[B A]
ans =
     1    -2    -2    -1     8     3
     3    -1     0    -2     3     1
     1     0    -3    -7     5     4
%This combined matrices B and A together, attaching A to the end of B.

[A x]
{_Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('horzcat')" style="font-weight:bold">horzcat</a>
Dimensions of matrices being concatenated are not consistent.}_ 

[A X]
ans =
    -2    -1     8     3    -2
     0    -2     3     1     3
    -3    -7     5     4     5
%This combined matrices A and X together, attaching X to the end of A.

[A ;y]
{_Error using <a href="matlab:matlab.internal.language.introspective.errorDocCallback('vertcat')" style="font-weight:bold">vertcat</a>
Dimensions of matrices being concatenated are not consistent.}_ 

[A ;x]
ans =
    -2    -1     8     3
     0    -2     3     1
    -3    -7     5     4
     1     2     3     4
%This added vector x to the bottom row of matrix A.


eye(5)
ans =
     1     0     0     0     0
     0     1     0     0     0
     0     0     1     0     0
     0     0     0     1     0
     0     0     0     0     1
%This returned an identity matrix with 5 rows and 5 columns

zeros(3,4)
ans =
     0     0     0     0
     0     0     0     0
     0     0     0     0
%This created a 3 by 4 empty matrix.

zeros(2)
ans =
     0     0
     0     0
%This created a 2 by 2 empty matrix.

ones(3,2)
ans =
     1     1
     1     1
     1     1
%This create a 3 by 2 matrix filled with ones.

ones(5)
ans =
     1     1     1     1     1
     1     1     1     1     1
     1     1     1     1     1
     1     1     1     1     1
     1     1     1     1     1
%This created a 5 by 5 matrix filled with ones.

diag([1 2 5 6 7])
ans =
     1     0     0     0     0
     0     2     0     0     0
     0     0     5     0     0
     0     0     0     6     0
     0     0     0     0     7
%This created a 5 by 5 matrix with the elements numbering down in a diagonal line.

diag([1 2 5 6 7],-1)
ans =
     0     0     0     0     0     0
     1     0     0     0     0     0
     0     2     0     0     0     0
     0     0     5     0     0     0
     0     0     0     6     0     0
     0     0     0     0     7     0
%This created the same matrix as the previous command, with the addition of an empty row on top and column on the right.

diag([1 2 5 6 7],2)
ans =
     0     0     1     0     0     0     0
     0     0     0     2     0     0     0
     0     0     0     0     5     0     0
     0     0     0     0     0     6     0
     0     0     0     0     0     0     7
     0     0     0     0     0     0     0
     0     0     0     0     0     0     0
%This replicated the previous matrix, with the addition of two empty columns to the left.

A,diag(A),diag(diag(A))
A =
    -2    -1     8     3
     0    -2     3     1
    -3    -7     5     4
ans =
    -2
    -2
     5
ans =
    -2     0     0
     0    -2     0
     0     0     5
%This returned matrix A, then a list of the diagonal elements of A, and finally a 3 by 3 matrix with A's diagonal elements going down.

magic(5)
ans =
    17    24     1     8    15
    23     5     7    14    16
     4     6    13    20    22
    10    12    19    21     3
    11    18    25     2     9
%This created a 5 by 5 matrix with random values

help magic
 <strong>magic</strong>  Magic square.
    <strong>magic</strong>(N) is an N-by-N matrix constructed from the integers
    1 through N^2 with equal row, column, and diagonal sums.
    Produces valid magic squares for all N > 0 except N = 2.

    <a href="matlab:doc magic">Reference page for magic</a>

%Gave a definition of the "magic" keyword and its applications.

hilb(5)
ans =
    1.0000    0.5000    0.3333    0.2500    0.2000
    0.5000    0.3333    0.2500    0.2000    0.1667
    0.3333    0.2500    0.2000    0.1667    0.1429
    0.2500    0.2000    0.1667    0.1429    0.1250
    0.2000    0.1667    0.1429    0.1250    0.1111
%Created a 5 by 5 matrix with each element being a fraction of 1 over the diagonal they're directly over.

help hilb
 <strong>hilb</strong>   Hilbert matrix.
    <strong>hilb</strong>(N) is the N by N matrix with elements 1/(i+j-1),
    which is a famous example of a badly conditioned matrix.
    See INVHILB for the exact inverse.
 
    <strong>hilb</strong>(N,CLASSNAME) produces a matrix of class CLASSNAME.
    CLASSNAME must be either 'single' or 'double' (the default).
 
    This is also a good example of efficient MATLAB programming
    style where conventional FOR or DO loops are replaced by
    vectorized statements. 
 
    See also <a href="matlab:help invhilb">invhilb</a>.

    <a href="matlab:doc hilb">Reference page for hilb</a>

%Gave a definition of the "hilb" command and its applications.


C=eye(3)
ans =
     1     0     0
     0     1     0
     0     0     1

D=diag([2 1 3],1)
ans =
     0     2     0     0
     0     0     1     0
     0     0     0     3
     0     0     0     0

E=ones(2,3)
ans =
     1     1     1
     1     1     1

V1=1:7
V1 =
     1     2     3     4     5     6     7
%This created a vector with elements starting at one and moving up in increments of one to seven.

V2=2:0.5:6.5
V2 =
    2.0000    2.5000    3.0000    3.5000    4.0000    4.5000    5.0000    5.5000    6.0000    6.5000
%This created a vector with elements starting at 2 and moving up in increments of 0.5 to 6.5.


V3=3:-1:-5
V3 =
     3     2     1     0    -1    -2    -3    -4    -5
%This created a vector with elements starting at 3 and moving down in decrements of 1 to -5.

V4=-5:1:1

V4 =

    -5    -4    -3    -2    -1     0     1

V5=10:-3:-2

V5 =

    10     7     4     1    -2

V6=5:-0.5:2

V6 =

    5.0000    4.5000    4.0000    3.5000    3.0000    2.5000    2.0000

V7 =

  Columns 1 through 9

         0    0.4000    0.8000    1.2000    1.6000    2.0000    2.4000    2.8000    3.2000

  Columns 10 through 11

    3.6000    4.0000


C;

C

C =

     1     0     0
     0     1     0
     0     0     1

R=434.1452

R =

  434.1452

%This created a variable R and assigned it the value 434.1452.

format long,R

R =

     4.341452000000000e+02

%This displayed R in scientific notation with 15 decimal places.
 
format short,R

R =

  434.1452

%This shortened R back to its initial 4 decimal places.

A,A+A

A =

    -2    -1     8     3
     0    -2     3     1
    -3    -7     5     4


ans =

    -4    -2    16     6
     0    -4     6     2
    -6   -14    10     8

A, 2*A

A =

    -2    -1     8     3
     0    -2     3     1
    -3    -7     5     4


ans =

    -4    -2    16     6
     0    -4     6     2
    -6   -14    10     8

A,B,A+B

A =

    -2    -1     8     3
     0    -2     3     1
    -3    -7     5     4


B =

     1    -2
     3    -1
     1     0

Matrix dimensions must agree.


B,E,B-2*E

B =

     1    -2
     3    -1
     1     0


E =

     1     1     1
     1     1     1

Matrix dimensions must agree.

x,X,x+X

x =

     1     2     3     4


X =

    -2
     3
     5


ans =

    -1     0     1     2
     4     5     6     7
     6     7     8     9

x,y,x+y

x =

     1     2     3     4


y =

     1
     4
     6
     3


ans =

     2     3     4     5
     5     6     7     8
     7     8     9    10
     4     5     6     7

%The errors that occur do so because the dimensions of the matrices must fit accordingly to add and multiply them and then combine.

transpose(A)

ans =

    -2     0    -3
    -1    -2    -7
     8     3     5
     3     1     4

A*D

ans =

     0    -4    -1    24
     0     0    -2     9
     0    -6    -7    15

A.*A

ans =

     4     1    64     9
     0     4     9     1
     9    49    25    16


A.*D
Matrix dimensions must agree.
 
A*A
Error using  * 
Inner matrix dimensions must agree.
 
%The errors occur because of mismatching matrix dimensions, such as one having less columns than the other.


G=[4 2 1
3 1 6
7 7 8]

G =

     4     2     1
     3     1     6
     7     7     8

G^2

ans =

    29    17    24
    57    49    57
   105    77   113

A^2
Error using  ^ 
One argument must be a square matrix and the other must be a scalar. Use POWER (.^) for
elementwise power.
 
%The command G^2 is equivalent to the command 'G,G*G'.
% The A^2 command resulted in an error, due to A not being a square matrix.
%The command 'A.*A' would calculate element-wise power 2 of A
 
 
rand(4)

ans =

    0.8147    0.6324    0.9575    0.9572
    0.9058    0.0975    0.9649    0.4854
    0.1270    0.2785    0.1576    0.8003
    0.9134    0.5469    0.9706    0.1419

rand(3,4)

ans =

    0.4218    0.9595    0.8491    0.7577
    0.9157    0.6557    0.9340    0.7431
    0.7922    0.0357    0.6787    0.3922

randi(100,2)

ans =

    66    71
    18     4

randi(10,2,4)

ans =

     3     1     7    10
     1     9     4     1

randi([10 40],2,4)

ans =

    23    33    15    23
    21    34    25    30

%rand(4) created a 4 by 4 matrix filled with random values.
%rand(4) created a 4 by 4 matrix filled with random values from 0 to 1.
%rand(3,4) creates a 3 by 4 matrix filled with random values from 0 to 1.
%randi(100,2) created a 2 by 2 matrix filled with random values from 0 to 100.
%randi(10,2,4) created a 2 by 4 matrix with random values from 0 to 10.
%randi([10 40],2,4) created a 2 by 4 matrix with random values from 10 to 40.
5*rand(3)

ans =

    3.5468    3.3985    0.5950
    3.7734    3.2755    2.4918
    1.3801    0.8131    4.7987

-3+5*rand(3)

ans =

   -1.2981    0.7563    0.4954
   -0.0737   -1.7245    1.4545
   -1.8809   -0.4702    1.7965

-3+5*rand(3)

ans =

   -0.2639   -1.7125    1.0714
   -2.3069    1.2036   -1.7824
   -2.2535   -1.7286    1.6463

randi([4 10],2,3)

ans =

     6     5     7
     5     8     6

randi([40 90],2,3)

ans =

    82    68    54
    69    86    78



diary on
%Exercise 3
type adds

function C=adds(A,B)
%This is the function that adds matrices A and B. It duplicates the MATLAB
%function A+B.
[m,n]=size(A);
[k,p]=size(B);
if m==k && n==p,
    for i=1:m
        for j=1:n
            C(i,j)=A(i,j)+B(i,j);
        end
    end
else
    disp('Error in using adds: matrices are not of the same size')
    C=[];
end
end

A=[1 2 3
4 5 6
7 8 9]

A =

     1     2     3
     4     5     6
     7     8     9

B=ones(2,3)

B =

     1     1     1
     1     1     1

C=adds(A,B)
Error in using adds: matrices are not of the same size

C =

     []

A=magic(3)

A =

     8     1     6
     3     5     7
     4     9     2

B=ones(3)

B =

     1     1     1
     1     1     1
     1     1     1

C=adds(A,B)

C =

     9     2     7
     4     6     8
     5    10     3

diary off

diary on
%Exercise 4

type sums

function A=sums(m,n)

A=[m,n];
if mod(m,1)~= 0 || mod(n,1) ~= 0
    A = [];
    disp('m and n must be integers')
    return,
end

for i=1:m
    for j=1:n
        A(i,j) = i+j;
    end
end


%(a)
m=6.5, n=4

m =

    6.5000


n =

     4

A=sums(m,n)
m and n must be integers

A =

     []

%(b)
m= 1, n=1.6

m =

     1


n =

    1.6000

>> A=sums(m,n)
m and n must be integers

A =

     []

%(c)
m=3, n=4

m =

     3


n =

     4

>> A=sums(m,n)

A =

     2     3     4     5
     3     4     5     6
     4     5     6     7


%(d)
m=3, n=3

m =

     3


n =

     3

>> A=sums(m,n)

A =

     2     3     4
     3     4     5
     4     5     6



function B=switches(A)
 
[m,n]=size(A);
 
for i=1:m
    for j=1:n
        B(j,i)=A(i,j);
    end
end


%(c)
A=sums(3,4)
A =

     2     3     4     5
     3     4     5     6
     4     5     6     7
B=switches(A=sums(3,4)

B =

     2     3     4
     3     4     5
     4     5     6
     5     6     7

%(d)
A=sums(3,3)
A =

     2     3     4
     3     4     5
     4     5     6
B=switches(A)

B =

     2     3     4
     3     4     5
     4     5     6


%The switches function does not actually alter the matrix output in part d, as both the columns and rows are the same within Matrix B and A’s square arrangement.

%The built-in ‘transpose’ function of MATLAB performs the same operation as that of the ‘switches’ function.

%(c transpose)
A=sums(3,4)

A =

     2     3     4     5
     3     4     5     6
     4     5     6     7

transpose(A)

ans =

     2     3     4
     3     4     5
     4     5     6
     5     6     7

%(d transpose)
A=sums(3,3)

A =

     2     3     4
     3     4     5
     4     5     6

transpose(A)

ans =

     2     3     4
     3     4     5
     4     5     6

diary off





