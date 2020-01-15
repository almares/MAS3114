%Exercise 1
diary Project3
diary on
format compact
type subspace

function[ ] = subspace(A,B)
m=size(A,1);
n=size(B,1);
%Subspace Calculation
if m~=n
    disp('Col A and Col B are subspaces of different spaces')
    return
else
    fprintf('Col A and Col B are subspaces of R^%i\n',m)
end
%Column Calculation
k = rank(A);
p = rank(B);
fprintf('The dim of Col A is k = %i\n',k);
fprintf('The dim of Col B is p = %i\n',p);
if k ~= p
    disp('k~=p, the dimensions of Col A and Col B are different')
else
    C = [A B];
    if rank(C) ~= k && rank(C) ~= p
        disp('k=p, the dimensions of Col A and Col B are the same, but Col    A ~= Col B')
    else
        disp('Col A = Col B')
    end
end
%All Real Calculation
if k == m
    fprintf('k=m (%i=%i) Col A is all R^%i\n',k,m,m);
else
    fprintf('k~=m (%i~=%i) Col A is not all R^%i\n',k,m,m);
end
if p == m
    fprintf('p=m (%i=%i) Col B is all R^%i\n',p,m,m);
else
    fprintf('p~=m (%i~=%i) Col B is not all R^%i\n',p,m,m);
end

A = [2 -4 -2 3; 6 -9 -5 8; 2 -7 -3 9; 4 -2 -2 -1; -6 3 3 4]
A =
     2    -4    -2     3
     6    -9    -5     8
     2    -7    -3     9
     4    -2    -2    -1
    -6     3     3     4
B = rref(A)
B =
    1.0000         0        -0.3333         0
         0      1.0000       0.3333         0
         0         0            0        1.0000
         0         0            0           0
         0         0            0           0
subspace(A,B)
Col A and Col B are subspaces of R^5
The dim of Col A is k = 3
The dim of Col B is p = 3
k=p, the dimensions of Col A and Col B are the same, but Col A ~= Col B
k~=m (3~=5) Col A is not all R^5
p~=m (3~=5) Col B is not all R^5
A = magic(4)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
B = eye(4)
B =
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
subspace(A,B)
Col A and Col B are subspaces of R^4
The dim of Col A is k = 3
The dim of Col B is p = 4
k~=p, the dimensions of Col A and Col B are different
k~=m (3~=4) Col A is not all R^4
p=m (4=4) Col B is all R^4
A = magic(4)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
B = eye(3)
B =
     1     0     0
     0     1     0
     0     0     1
subspace(A,B)
Col A and Col B are subspaces of different spaces
A = magic(5)
A =
    17    24     1     8    15
    23     5     7    14    16
     4     6    13    20    22
    10    12    19    21     3
    11    18    25     2     9
B = eye(5)
B =
     1     0     0     0     0
     0     1     0     0     0
     0     0     1     0     0
     0     0     0     1     0
     0     0     0     0     1
subspace(A,B)
Col A and Col B are subspaces of R^5
The dim of Col A is k = 5
The dim of Col B is p = 5
Col A = Col B
k=m (5=5) Col A is all R^5
p=m (5=5) Col B is all R^5
% The elementary row operations might change the column space since the row reduced form of a
% matrix has different columns than a non-reduced form of that matrix, and therefore a different 
% column space
diary off

% Exercise 2
diary on
format compact
type shrink

function B=shrink(A)
[~,pivot]=rref(A);
B=A(:,pivot);
end

A=magic(4)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
[~,pivot]=rref(A)
pivot =
     1     2     3
% This command gives which columns of the matrix A are pivot columns
B=A(:,pivot)
B =
    16     2     3
     5    11    10
     9     7     6
     4    14    15
% This command gives columns from A which are pivot columns. We can see that matrix B is the columns 1,2, and 3 of A.

type basis

function B=basis(A)
m=size(A,1);
A=shrink(A);
sprintf('A basis for Col A is \n')
B=A
x=rank(B);
    if x==m
        sprintf('A basis for R^%i is \n',m)
    else D=eye(m);
        D=[B D];
        D=shrink(D);
        y=rank(D);
        if y==m
            sprintf('A basis for R^%i is \n',m)
            B=D;
        else disp('What? It is not a basis!?')
        end
    end
end
A=[1 0;0 0;0 0;0 1]
A =
     1     0
     0     0
     0     0
     0     1
B=basis(A)
ans =
A basis for Col A is 

B =
     1     0
     0     0
     0     0
     0     1
ans =
A basis for R^4 is 

B =
     1     0     0     0
     0     0     1     0
     0     0     0     1
     0     1     0     0
A=[2 0;4 0;1 0;0 0]
A =
     2     0
     4     0
     1     0
     0     0
B=basis(A)
ans =
A basis for Col A is 

B =
     2
     4
     1
     0
ans =
A basis for R^4 is 

B =
     2     1     0     0
     4     0     1     0
     1     0     0     0
     0     0     0     1
A=magic(3)
A =
     8     1     6
     3     5     7
     4     9     2
B=basis(A)
ans =
A basis for Col A is 

B =
     8     1     6
     3     5     7
     4     9     2
ans =
A basis for R^3 is 

B =
     8     1     6
     3     5     7
     4     9     2
A=magic(6)
A =
    35     1     6    26    19    24
     3    32     7    21    23    25
    31     9     2    22    27    20
     8    28    33    17    10    15
    30     5    34    12    14    16
     4    36    29    13    18    11
B=basis(A)
ans =
A basis for Col A is 

B =
    35     1     6    26    19
     3    32     7    21    23
    31     9     2    22    27
     8    28    33    17    10
    30     5    34    12    14
     4    36    29    13    18
ans =
A basis for R^6 is 

B =
    35     1     6    26    19     1
     3    32     7    21    23     0
    31     9     2    22    27     0
     8    28    33    17    10     0
    30     5    34    12    14     0
     4    36    29    13    18     0
diary off









%Exercise 3
format compact
diary on
type closetozeroroundoff

function B = closetozeroroundoff(A)
[m,n]=size(A);
for i=1: m
    for j=1: n
        if abs(A(i,j))<10^(-7)
            A(i,j)=0;
        end
    end
end
B=A;

 type polyspace

function B = polyspace(P,Q,r)
format rat,
u = sym2poly(P(1));
n = length(u);

C = zeros(n);
for i=1: n
    C(:,i)=transpose(sym2poly(P(i)));
end
B=closetozeroroundoff(C);

x = rank(B);

if x~=n
    sprintf('The polynomials in P do not form a basis for P%d',n-1)
    fprintf('The reduced echelon form of B is %\n')
    A = rref(B);
    return
else
    sprintf('The polynomials in P do form a basis for P%d',n-1)
    v = sym2poly(Q(1));
    v = transpose(v);
    y = [B v];
    y = rref(y);
    y = y(:,n+1);
    fprintf('The coordinates of the polynomial Q with respect to the basis P are %\n')
    y = closetozeroroundoff(y)
    
    fprintf('R, given vector r is:')
    R = r/v
end
end

 syms x;
 
 %a)
P = [x^3+(3*x^2), 10^(-8)*x^3+x, 10^(-8)*x^3+4*x^2+x,x^3+x];
Q = [10^(-8)*x^3 + x^2 + 6*x];
r = [2;-3;1;0];
B = polyspace(P,Q,r)
ans =
    'The polynomials in P do not form a basis for P3'
The reduced echelon form of B is B =
       1              0              0              1       
       3              0              4              0       
       0              1              1              1       
       0              0              0              0       

 %b)
 P = [x^3-1,10^(-8)*x^3+2*x^2, 10^(-8)*x^3+x,x^3+x];
 B = polyspace(P,Q,r)
ans =
    'The polynomials in P do form a basis for P3'
The coordinates of the polynomial Q with respect to the basis P are y =
       0       
       1/2     
       6       
       0       
R, given vector r is:R =
       0              0              1/3            0       
       0              0             -1/2            0       
       0              0              1/6            0       
       0              0              0              0       
B =
       1              0              0              1       
       0              2              0              0       
       0              0              1              1       
      -1              0              0              0       

 %c) 
 P = [x^4+x^3+x^2+1, 10^(-8)*x^4+x^3+x^2+x+1, 10^(-8)*x^4+x^2+x+1, 10^(-8)*x^4+x+1, 10^(-8)*x^4+1];
Q = [x^4-1];
r = diag(magic(5));
B = polyspace(P,Q,r)
ans =
    'The polynomials in P do form a basis for P4'
The coordinates of the polynomial Q with respect to the basis P are y =
       1       
      -1       
       0       
       1       
      -2       
R, given vector r is:R =
      17              0              0              0              0       
       5              0              0              0              0       
      13              0              0              0              0       
      21              0              0              0              0       
       9              0              0              0              0       
B =
       1              0              0              0              0       
       1              1              0              0              0       
       1              1              1              0              0       
       0              1              1              1              0       
       1              1              1              1              1        
 diary off









%Exercise 4

type reimsum

function [T,I] = reimsum( P,a,b,n )
%REIMSUM This function take polynomial P, scalars a 
%and b, and vector n and uses Riemann sums to approximate
%the definite integral on the interval [a,b] by
%subintervals of equal length h=(b-a)/n(j)

format long

n=transpose(n);
c = zeros(length(n),1);
d = zeros(length(n),1);
f = zeros(length(n),1);
h=(b-a)./n;
s2p = sym2poly(P);

%c (leftmost)

for i=1:length(n)
    reimann_sum=0;
    for j=1:n(i)
        left_point = a+(j-1)*h(i);
        reimann_sum = reimann_sum+polyval(s2p, left_point)*h(i);
    end
    c(i) = reimann_sum;
end

%d (center)
for i=1:length(n)
    reimann_sum = 0;
    for j=1:n(i)
        center_point = a+(j-1/2)*h(i);
        reimann_sum = reimann_sum + polyval(s2p, center_point)*h(i);
    end
    d(i) = reimann_sum;
end
d = closetozeroroundoff(d);

%f (rightmost)
for i=1:length(n)
    reimann_sum = 0;
    for j=1:n(i)
        right_point = a+(j)*h(i);
        reimann_sum = reimann_sum + polyval(s2p, right_point)*h(i);
    end
    f(i) = reimann_sum;
end

A = [n c d f];
T = array2table(A, 'VariableNames', {'n', 'Left', 'Middle', 'Right'});

I = double(int(P, a, b));

%(a) 
P=2*x^4+4*x^2-1
 
P =
 
2*x^4 + 4*x^2 - 1
 
a=-1; b = 1;
n=[1:10]

n =

     1     2     3     4     5     6     7     8     9    10

[T,I] = reimsum(P,a,b,n)

T = 

Left        	 Middle      		Right     
______________________________________________________________

     1                  10                   -2                  10
     2                   4                 0.25                   4
     3    2.62551440329218    0.897119341563786    2.62551440329218
     4               2.125             1.140625               2.125
     5             1.88992              1.25632             1.88992
     6    1.76131687242798     1.31995884773663    1.76131687242798
     7    1.68346522282382     1.35860058309038    1.68346522282382
     8           1.6328125         1.3837890625           1.6328125
     9     1.5980287557791     1.40110755474267     1.5980287557791
    10             1.57312              1.41352             1.57312


I =

   1.466666666666667

n=[1,5,10,100,1000,10000]

n =

  Columns 1 through 5

           1           5          10         100        1000

  Column 6

       10000

[T,I] = reimsum(P,a,b,n)




T = 

Left        	 Middle      		Right     
______________________________________________________________

        1                  10                  -2                  10
        5             1.88992             1.25632             1.88992
       10             1.57312             1.41352             1.57312
      100         1.467733312         1.466133352         1.467733312
     1000     1.4666773333312     1.4666613333352     1.4666773333312
    10000    1.46666677333334    1.46666661333333    1.46666677333334


I =

   1.466666666666667

a=-10; b = 10;
[T,I] = reimsum(P,a,b,n)

T = 

 Left        	 Middle      		Right     
______________________________________________________________

        1              407980                 -20              407980
        5              103852               72172              103852
       10               88012               79972               88012
      100          82700.5312          82619.7352          82700.5312
     1000    82647.2053331201      82646.39733352    82647.2053331201
    10000    82646.6720533329    82646.6639733331    82646.6720533329


I =

     8.264666666666667e+04

n=[1:10]

n =

     1     2     3     4     5     6     7     8     9    10

[T,I] = reimsum(P,a,b,n)
 


T = 

Left        	 Middle      		Right     
______________________________________________________________

     1              407980                 -20              407980
     2              203980               26980              203980
     3    139864.773662551    55025.2674897119    139864.773662551
     4              115480             66542.5              115480
     5              103852               72172              103852
     6    97445.0205761317    75309.2181069959    97445.0205761317
     7    93551.0120783007    77227.8134110787    93551.0120783007
     8            91011.25         78483.90625            91011.25
     9    89264.3570593914    79350.0147335264    89264.3570593914
    10               88012               79972               88012


I =

     8.264666666666667e+04

%(b)
P=x^3-2*x
 
P =
 
x^3 - 2*x
 
a=-1; b = 1;
n=[1:10]

n =

     1     2     3     4     5     6     7     8     9    10

[T,I] = reimsum(P,a,b,n)

T = 

Left         Middle      		Right     
______________________________________________________________

     1                    2    0                         -2
     2                    1    0                         -1
     3    0.666666666666667    0         -0.666666666666667
     4                  0.5    0                       -0.5
     5                  0.4    0                       -0.4
     6    0.333333333333333    0         -0.333333333333333
     7    0.285714285714286    0         -0.285714285714285
     8                 0.25    0                      -0.25
     9    0.222222222222222    0         -0.222222222222222
    10                  0.2    0                       -0.2


I =

     0

n=[1,5,10,100,1000,10000]

n =

  Columns 1 through 5

           1           5          10         100        1000

  Column 6

       10000

[T,I] = reimsum(P,a,b,n)

T = 

Left         Middle      		Right     
______________________________________________________________

        1                       2    0                            -2
        5                     0.4    0                          -0.4
       10                     0.2    0                          -0.2
      100      0.0200000000000001    0           -0.0200000000000001
     1000     0.00199999999999985    0          -0.00199999999999995
    10000    0.000199999999999973    0         -0.000200000000000095


I =

     0

a=-10; b = 10;
[T,I] = reimsum(P,a,b,n)

T = 

Left         Middle      		Right     
______________________________________________________________

        1               -19600    0                    19600
        5                -3920    0                     3920
       10                -1960    0                     1960
      100    -195.999999999999    0         196.000000000001
     1000    -19.6000000000001    0         19.6000000000003
    10000    -1.96000000000002    0         1.95999999999917


I =

     0

n=[1:10]

n =

     1     2     3     4     5     6     7     8     9    10

[T,I] = reimsum(P,a,b,n)

T = 

Left       	 Middle      		Right     
______________________________________________________________

     1               -19600    0                    19600
     2                -9800    0                     9800
     3    -6533.33333333333    0         6533.33333333333
     4                -4900    0                     4900
     5                -3920    0                     3920
     6    -3266.66666666667    0         3266.66666666667
     7                -2800    0                     2800
     8                -2450    0                     2450
     9    -2177.77777777778    0         2177.77777777778
    10                -1960    0                     1960


I =

     0

%The middle point seems to be the best approximation, especially
%when I=0
diary off










%Exercise 5
type polint

function B = polint( P )
%POLINT accepts a polynomial P and outputs an antiderivative 
%of the polynomial B
syms x;

%the new coefficient is equal to the old coefficient divided
%by its exponent + 1
for n=1:length(P)
    coeff(n) = (P(n)/(length(P)-n+1));
end

%Assign 5 to arbitrary constant 
coeff(length(P)+1) = 5;
B = poly2sym(coeff,x);
end

%(a)
P=[6,5,4,3,2,6]
P =
     6     5     4     3     2     6
B=polint(P)
B =
x^6 + x^5 + x^4 + x^3 + x^2 + 6*x + 5
int(poly2sym(P))+5
ans =
x^6 + x^5 + x^4 + x^3 + x^2 + 6*x + 5'

%(b)
P=[1,0,-1,0,3,0,1]
P =
     1     0    -1     0     3     0     1
B=polint(P)
B =
x^7/7 - x^5/5 + x^3 + x + 5
int(poly2sym(P))+5
ans =
x^7/7 - x^5/5 + x^3 + x + 5
diary off








%Exercise6 
diary on
type markov
 
function q = markov(P,x0)
k = 0; 
xk = x0; 
n = size(P,1);
s = (sum(P,1));
S = ones(1,n);
    if s == S
        Q = null(P-eye(n),'r');
        c = sum(Q);
        q = (1./c).*Q;
    while norm(xk-q) > 10^(-7)
        xk = P*xk;
        k = k + 1;
    end
        k
        xk
    else
         disp('P is not a stochastic matrix')
         q = [];
    end
end

%a 
P = [.6,.3;.5,.7];
x0 = [.4;.6];
q = markov(P,x0)
P is not a stochastic matrix
q =
     []
 
%b
P = [.5,.3;.5,.7];
 
q = markov(P,x0)
k =
     8
xk =
    0.3750
    0.6250
q =
    0.3750
    0.6250
% xk and q are the same vector 
 
%c
P = [.9,.2;.1,.8];
x0 = [.12;.88];
q = markov(P,x0)
k =
    45
xk =
    0.6667
    0.3333
q =
    0.6667
    0.3333
 
%d
x0 = [.14;.86];
q = markov(P,x0)
k =
    45
xk =
    0.6667
    0.3333
q =
    0.6667
    0.3333
x0 = [.86;.14];
q = markov(P,x0)
k =
    42
xk =
    0.6667
    0.3333
q =
    0.6667
    0.3333
 
% q in part (d) and (c) is the same 
%the output is not affected by different x0 vector elements or order as long as those element add up to one
%the elements can contribute to a change in the number of k (interation) 
 
%e
P = [.9,.01,.09;.01,.9,.01;.09,.09,.9];
x0 = [.5;.3;.2];
q = markov(P,x0)
k =
   128
xk =
    0.4354
    0.0909
    0.4737
q =
    0.4354
    0.0909
    0.4737
diary off

