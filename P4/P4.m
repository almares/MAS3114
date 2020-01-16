diary Project4
format compact

%Exercise1
type eigen

function [] = eigen(A)
%EIGEN will find all the eigenvalues of a nxn matrix,  find 
%a rational basis for each eigenvalue, decide if the matrix
%is diagonalizable, and if it is, output the invertible 
%matrix P and diagonal matrix D

%1: 
L=closetozeroroundoff(transpose(eig(A)))

%2: 
M=unique(L)

%3:
num_eig=size(L);
disp 'The sum of multiplicities of the eigenvalues is'
Q=num_eig(1,2)

%4: 
size_unique=size(M);
num_unique=size_unique(1,2);
for i=1:num_unique
    W=null(A-M(1,i)*eye(size(A)),'r')
end

%5:
dim_num=null(A-M(1,1)*eye(size(A)),'r');
for j=2:num_unique
    dim_num=[dim_num null(A-M(1,j)*eye(size(A)),'r')];
end 
size_dim_num=size(dim_num);
disp 'The total sum of the dimensions of the eigenspace is'
N=size_dim_num(1,2)

if N<Q
    disp 'No, Matrix A is not diagonalizable since N<Q'
else
    disp 'Yes, Matrix A is diagonalizable since N=Q'
    
%6: 
    P=dim_num
    D=diag(L)
    F=closetozeroroundoff(A*P-P*D)
    if F==0
        disp 'Great! I got a diagonalization!'
    else
        disp 'Oops! I got a bug in my code'
    end
end
end



type closetozeroroundoff

function B = closetozeroroundoff(A)
%CLOSETOZEROROUNDOFF 
[m,n]=size(A);
for i=1:m
    for j=1:n
        if abs(A(i,j))<10^(-7)
            A(i,j) = 0;
        end
    end
end
B=A;

end

type jord

function J = jord(n,r)
A = ones(n);
J = tril(triu(A),1);
for i = 1: n
    J(i, i) = r;
end

%(a)
A=[2 2; 0 2];
eigen(A)

L =
     2     2
M =
     2

The sum of multiplicities of the eigenvalues is
Q =
     2
W =
     1
     0

The total sum of the dimensions of the eigenspace is
N =
     1

No, Matrix A is not diagonalizable since N<Q

%(b)
A=[4 0 0 0; 1 3 0 0; 0 -1 3 0; 0 -1 5 4];
eigen(A)

L =
     4     3     3     4

M =
     3     4
The sum of multiplicities of the eigenvalues is
Q =
     4
W =
         0
         0
   -0.2000
    1.0000
W =

     0
     0
     0
     1

The total sum of the dimensions of the eigenspace is
N =
     2

No, Matrix A is not diagonalizable since N<Q
%(c)
A=jord(5,3);
eigen(A)
L =
     3     3     3     3     3
M =
     3
The sum of multiplicities of the eigenvalues is
Q =
     5
W =
     1
     0
     0
     0
     0

The total sum of the dimensions of the eigenspace is
N =
     1
No, Matrix A is not diagonalizable since N<Q
%(d)
A=diag([3, 3, 3, 2, 2, 1]);
eigen(A)

L =
     1     2     2     3     3     3
M =
     1     2     3
The sum of multiplicities of the eigenvalues is
Q =
     6
W =
     0
     0
     0
     0
     0
     1
W =

     0     0
     0     0
     0     0
     1     0
     0     1
     0     0
W =

     1     0     0
     0     1     0
     0     0     1
     0     0     0
     0     0     0
     0     0     0
The total sum of the dimensions of the eigenspace is
N =
     6
Yes, Matrix A is diagonalizable since N=Q
P =
     0     0     0     1     0     0
     0     0     0     0     1     0
     0     0     0     0     0     1
     0     1     0     0     0     0
     0     0     1     0     0     0
     1     0     0     0     0     0
D =
     1     0     0     0     0     0
     0     2     0     0     0     0
     0     0     2     0     0     0
     0     0     0     3     0     0
     0     0     0     0     3     0
     0     0     0     0     0     3
F =
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0
     0     0     0     0     0     0

Great! I got a diagonalization!

%(e)
A=magic(4);
eigen(A)

L =
   -8.9443         0    8.9443   34.0000
M =
   -8.9443         0    8.9443   34.0000

The sum of multiplicities of the eigenvalues is
Q =
     4
W =
   -0.4570
   -0.0287
   -0.5143
    1.0000
W =
    -1
    -3
     3
     1
W =
   -2.1882
    1.1254
    0.0627
    1.0000
W =
    1.0000
    1.0000
    1.0000
    1.0000

The total sum of the dimensions of the eigenspace is
N =
     4

Yes, Matrix A is diagonalizable since N=Q

P =
   -0.4570   -1.0000   -2.1882    1.0000
   -0.0287   -3.0000    1.1254    1.0000
   -0.5143    3.0000    0.0627    1.0000
    1.0000    1.0000    1.0000    1.0000
D =
   -8.9443         0         0         0
         0         0         0         0
         0         0    8.9443         0
         0         0         0   34.0000
F =
     0     0     0     0
     0     0     0     0
     0     0     0     0
     0     0     0     0

Great! I got a diagonalization!

%(f)
A=ones(5);
eigen(A)

L =
     0     0     0     0     5
M =
     0     5

The sum of multiplicities of the eigenvalues is
Q =
     5
W =
    -1    -1    -1    -1
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
W =
     1
     1
     1
     1
     1

The total sum of the dimensions of the eigenspace is
N =
     5

Yes, Matrix A is diagonalizable since N=Q
P =
    -1    -1    -1    -1     1
     1     0     0     0     1
     0     1     0     0     1
     0     0     1     0     1
     0     0     0     1     1
D =
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     5
F =
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0

Great! I got a diagonalization!

%(g) 
A=magic(5);
eigen(A)

L =
  -21.2768  -13.1263   13.1263   21.2768   65.0000
M =
  -21.2768  -13.1263   13.1263   21.2768   65.0000

The sum of multiplicities of the eigenvalues is
Q =
     5
W =
   Empty matrix: 5-by-0
W =
   -2.4172
    2.2511
   -1.4952
    0.6613
    1.0000
W =
   -0.4137
   -0.2736
    0.6186
   -0.9313
    1.0000
W =
   Empty matrix: 5-by-0
W =
   Empty matrix: 5-by-0

The total sum of the dimensions of the eigenspace is
N =
     2

No, Matrix A is not diagonalizable since N<Q

%The eigen function isn't correct for g because the built in 
%matlab function for diagonalization shows that it is 
%diagonizable

eigen_1(A)

L =
  -21.2768  -13.1263   13.1263   21.2768   65.0000
M =
  -21.2768  -13.1263   13.1263   21.2768   65.0000

The sum of multiplicities of the eigenvalues is
Q =
     5
W =
   -0.0976
   -0.3525
   -0.5501
    0.3223
    0.6780


W =
   -0.6330
    0.5895
   -0.3915
    0.1732
    0.2619
W =
   -0.2619
   -0.1732
    0.3915
   -0.5895
    0.6330
W =
    0.6780
    0.3223
   -0.5501
   -0.3525
   -0.0976
W =
   -0.4472
   -0.4472
   -0.4472
   -0.4472
   -0.4472

The total sum of the dimensions of the eigenspace is
N =
     5

Yes, Matrix A is diagonalizable since N=Q
P =
   -0.0976   -0.6330   -0.2619    0.6780   -0.4472
   -0.3525    0.5895   -0.1732    0.3223   -0.4472
   -0.5501   -0.3915    0.3915   -0.5501   -0.4472
    0.3223    0.1732   -0.5895   -0.3525   -0.4472
    0.6780    0.2619    0.6330   -0.0976   -0.4472
D =
  -21.2768         0         0         0         0
         0  -13.1263         0         0         0
         0         0   13.1263         0         0
         0         0         0   21.2768         0
         0         0         0         0   65.0000
F =
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0
     0     0     0     0     0

Great! I got a diagonalization!
%I changed the rational basis in part 3 to use null() instead 
%of null( , 'r'). The code is now correct because it shows that
%g is diagonalizable.




























diary on
format compact
%Exercise 2
type diagonal
function L = diagonal(A)
n = size(A,1)
[P,D] = eig(A);
k = rank(P);
fprintf('The number of linearly independent columns in P is k = %g\n',k)
if k == n
    disp('A is diagonalizable')
    disp('A basis for R^n is')
    disp(P)
else
    disp('A is not diagonalizable')
    disp('A does not have enough linearly independent eigenvectors to create a basis for R^n')
end
L = transpose(diag(D));
%a
 A = [2 2;0 2];
L = diagonal(A)
n =
     2
The number of linearly independent columns in P is k = 1
A is not diagonalizable
A does not have enough linearly independent eigenvectors to create a basis for R^n
L =
     2     2
%b
A = [4 0 0 0;1 3 0 0; 0 -1 3 0; 0 -1 5 4];
L = diagonal(A)
n =
     4
The number of linearly independent columns in P is k = 2
A is not diagonalizable
A does not have enough linearly independent eigenvectors to create a basis for R^n
L =
     4     3     3     4
%c
A = jord(5,3);
L = diagonal(A)
n =
     5
The number of linearly independent columns in P is k = 1
A is not diagonalizable
A does not have enough linearly independent eigenvectors to create a basis for R^n
L =
     3     3     3     3     3
%d
A = diag([3,3,3,2,2,1]);
L = diagonal(A)
n =
     6
The number of linearly independent columns in P is k = 6
A is diagonalizable
A basis for R^n is
     0     0     0     0     0     1
     0     0     0     1     0     0
     0     0     0     0     1     0
     0     1     0     0     0     0
     0     0     1     0     0     0
     1     0     0     0     0     0
L =
     1     2     2     3     3     3
%e
A = magic(4);
L = diagonal(A)
n =
     4
The number of linearly independent columns in P is k = 4
A is diagonalizable
A basis for R^n is
   -0.5000   -0.8236    0.3764   -0.2236
   -0.5000    0.4236    0.0236   -0.6708
   -0.5000    0.0236    0.4236    0.6708
   -0.5000    0.3764   -0.8236    0.2236
L =
   34.0000    8.9443   -8.9443    0.0000
%f
A = ones(5);
L = diagonal(A)
n =
     5
The number of linearly independent columns in P is k = 5
A is diagonalizable
A basis for R^n is
    0.8333   -0.1667   -0.1667    0.2236    0.4472
   -0.1667    0.8333   -0.1667    0.2236    0.4472
   -0.1667   -0.1667    0.8333    0.2236    0.4472
   -0.5000   -0.5000   -0.5000    0.2236    0.4472
         0         0         0   -0.8944    0.4472
L =
     0     0     0     0     5
%g
A = magic(5);
L = diagonal(A)
n =
     5
The number of linearly independent columns in P is k = 5
A is diagonalizable
A basis for R^n is
   -0.4472    0.0976   -0.6330    0.6780   -0.2619
   -0.4472    0.3525    0.5895    0.3223   -0.1732
   -0.4472    0.5501   -0.3915   -0.5501    0.3915
   -0.4472   -0.3223    0.1732   -0.3525   -0.5895
   -0.4472   -0.6780    0.2619   -0.0976    0.6330
L =
   65.0000  -21.2768  -13.1263   21.2768   13.1263














diary Project4
diary on
format compact
%Exercise 3
type shrink
function B = shrink(A)
format compact,
[~,pivot]=rref(A);
B=A(: ,pivot);
end
type proj
function p = proj(A,b)
format compact
A=shrink(A);
b=transpose(b);
At=A';
m=length(At(1,:));
if length(b)~=m
    disp('No solution: dimensions of A and B disagree')
    disp('p=empty_vector')
    return
end
if rank([A b])==rank(A)
    disp('b is in ColA');
    return
end
if A*transpose(A)*b == zeros(m,1)
    p=zeros(m,1);
    z=b;
    disp('b is orthogonal to ColA')
    return
end
xthat=(A'*A)\(A'*b);
p=A*xthat;
z=b-p;
if abs(p'*z)<=10^-7
    disp('Yes! p and z are orthogonal')
else
    disp('Oops! Is there a bug in my code?')
end
end

%a)
A = magic(6); A = A(:,1:4),b=(1:6)
A =
    35     1     6    26
     3    32     7    21
    31     9     2    22
     8    28    33    17
    30     5    34    12
     4    36    29    13
b =
     1     2     3     4     5     6
proj(A,b)
Yes! p and z are orthogonal
ans =
    0.9492
    2.1599
    2.9492
    3.9180
    5.1287
    5.9180

%b)
A = magic(6), E = eye(6); b = E(6,:)
A =
    35     1     6    26    19    24
     3    32     7    21    23    25
    31     9     2    22    27    20
     8    28    33    17    10    15
    30     5    34    12    14    16
     4    36    29    13    18    11
b =
     0     0     0     0     0     1
proj(A,b)
Yes! p and z are orthogonal
ans =
   -0.2500
   -0.0000
    0.2500
    0.2500
   -0.0000
    0.7500
 
%c)
A = magic(4), b = (1:5)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
b =
     1     2     3     4     5
proj(A,b)
No solution: dimensions of A and B disagree
p=empty_vector

%d)
A = magic(5), b = rand(1,5)
A =
    17    24     1     8    15
    23     5     7    14    16
     4     6    13    20    22
    10    12    19    21     3
    11    18    25     2     9
b =
    0.8147    0.9058    0.1270    0.9134    0.6324
proj(A,b)
b is in ColA
%e)
A = ones(6); A (:) = 1:36, b = [1,0,1,0,1,0]
A =
     1     7    13    19    25    31
     2     8    14    20    26    32
     3     9    15    21    27    33
     4    10    16    22    28    34
     5    11    17    23    29    35
     6    12    18    24    30    36
b =
     1     0     1     0     1     0
proj(A,b)
Yes! p and z are orthogonal
ans =
    0.7143
    0.6286
    0.5429
    0.4571
    0.3714
    0.2857
%f)
A = ones(6);A(:)=1:36; A = null(A,'r'), b = ones(1,6)
A =
     1     2     3     4
    -2    -3    -4    -5
     1     0     0     0
     0     1     0     0
     0     0     1     0
     0     0     0     1
b =
     1     1     1     1     1     1
proj(A,b)
b is orthogonal to ColA
ans =
     0
     0
     0
     0
     0
     0

diary off




diary on
format compact
%Exercise 4
type shrink

function B = shrink(A)
format compact
[~,pivot] = rref(A);
B = A(:,pivot);
type solvemore

function X = solvemore(A,b)
format long,
A=shrink(A);
[m,n]=size(A);
[p,~]=size(b);
B = closetozeroroundoff(A\A-eye(n));
if p == rank(A)
    disp('The equation is consistent - look for the exact solution')
    x1 = A\b;
    if B ~= zeros(n)
        disp('A does not have orthonormal columns')
        X = x1;
    elseif m ~= n
        disp('A has orthonormal columns but is not orthogonal')
        X = x1;
    else
        disp('A is orthogonal')
        x2 = inv(A)*b;
        N = norm(x1-x2);
        X = [x1,x2];
        fprintf('The norm of the difference between two solutions is N = %g\n',N)
    end
else
    disp('The system is inconsistent - look for the least-squares solution')
    x3 = A\(A*(A\b));
    n1 = norm(b-A*x3);
    fprintf('The solution of the normal equations is\n')
    disp(x3)
    fprintf('The least-squares error of the approximation is n1 = %g\n',n1)
    if B==zeros(n)
        disp('A has orthonormal columns: an orthonormal basis for Col A is U=A')
        U = A;
    else
        U = orth(A);
        disp('An orthonormal basis for Col A is U =')
        disp(U)
    end
    b1 = U/U*b;
    disp('The projection of b onto Col A is')
    disp(b1)
    x4 = A\(b1);
    disp('The least-squares solution by using the projection onto Col A is x4 =')
    disp(x4)
    n2 = norm(b-b1);
    fprintf('The least-squares error of this approximation is n2 = %g\n',n2)
    n3 = norm(x3-x4);
    fprintf('The least-squares error of this approximation is n3 = %g\n',n3)
    x = rand(n,1);
    n4 = norm(b-A*x);
    fprintf('An error of approximation of b by Ax for a random vector x in R^n is %g\n',n4)
    X = [x3, x4];
End

%A
A=magic(4); b=A(:,1),A=orth(A)
b =
    16
     5
     9
     4
A =
  -0.500000000000000   0.670820393249937   0.500000000000000
  -0.500000000000000  -0.223606797749979  -0.500000000000000
  -0.500000000000000   0.223606797749979  -0.500000000000000
  -0.500000000000000  -0.670820393249937   0.500000000000000

X = solvemore(A,b)
The system is inconsistent - look for the least-squares solution
The solution of the normal equations is
 -17.000000000000000
   8.944271909999165
   3.000000000000002
The least-squares error of the approximation is n1 = 1.4349e-14
A has orthonormal columns: an orthonormal basis for Col A is U=A
The projection of b onto Col A is
  16.000000000000004
   5.000000000000002
   9.000000000000000
   4.000000000000001
The least-squares solution by using the projection onto Col A is x4 =
 -16.999999999999996
   8.944271909999161
   3.000000000000001
The least-squares error of this approximation is n2 = 4.07014e-15
The least-squares error of this approximation is n3 = 5.1022e-15
An error of approximation of b by Ax for a random vector x in R^n is 19.96
X =
 -17.000000000000000 -16.999999999999996
   8.944271909999165   8.944271909999161
   3.000000000000002   3.000000000000001

%B
A=magic(5);A=orth(A),b=rand(5,1)
A =
  Columns 1 through 3
  -0.447213595499958  -0.545634873129948   0.511667273601714
  -0.447213595499958  -0.449758363151205  -0.195439507584838
  -0.447213595499958  -0.000000000000024  -0.632455532033676
  -0.447213595499958   0.449758363151189  -0.195439507584872
  -0.447213595499958   0.545634873129987   0.511667273601672
  Columns 4 through 5
   0.195439507584854  -0.449758363151198
  -0.511667273601691   0.545634873129969
   0.632455532033676  -0.000000000000002
  -0.511667273601694  -0.545634873129966
   0.195439507584856   0.449758363151196
b =
   0.278498218867048
   0.546881519204984
   0.957506835434298
   0.964888535199277
   0.157613081677548

X = solvemore(A,b)
The equation is consistent - look for the exact solution
A is orthogonal
The norm of the difference between two solutions is N = 3.79044e-16
X =
  -1.299329098944367  -1.299329098944367
   0.122043004805592   0.122043004805592
  -0.677896209908240  -0.677896209908240
  -0.082709389188605  -0.082709389188605
  -0.282448306571140  -0.282448306571140

%C
A=magic(4),b=ones(4,1)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
b =
     1
     1
     1
     1

X = solvemore(A,b)
The system is inconsistent - look for the least-squares solution
The solution of the normal equations is
   0.058823529411765
   0.117647058823529
  -0.058823529411765
The least-squares error of the approximation is n1 = 7.19507e-16
A has orthonormal columns: an orthonormal basis for Col A is U=A
The projection of b onto Col A is
   1.000000000000000
   1.000000000000000
   1.000000000000000
   1.000000000000000
The least-squares solution by using the projection onto Col A is x4 =
   0.058823529411765
   0.117647058823529
  -0.058823529411765
The least-squares error of this approximation is n2 = 4.15407e-16
The least-squares error of this approximation is n3 = 3.58211e-16
An error of approximation of b by Ax for a random vector x in R^n is 39.3267
X =
   0.058823529411765   0.058823529411765
   0.117647058823529   0.117647058823529
  -0.058823529411765  -0.058823529411765

%D
A=magic(4),b=rand(4,1)
A =
    16     2     3    13
     5    11    10     8
     9     7     6    12
     4    14    15     1
b =
   0.800280468888800
   0.141886338627215
   0.421761282626275
   0.915735525189067

X = solvemore(A,b)
The system is inconsistent - look for the least-squares solution
The solution of the normal equations is
   0.032693337120072
  -0.221874291154678
   0.256229883897294
The least-squares error of the approximation is n1 = 0.213562
A has orthonormal columns: an orthonormal basis for Col A is U=A
The projection of b onto Col A is
   0.800280468888800
   0.460246301393030
   0.421761282626275
   0.915735525189066
The least-squares solution by using the projection onto Col A is x4 =
   0.033395601743820
  -0.136900271681156
   0.179917128116664
The least-squares error of this approximation is n2 = 0.31836
The least-squares error of this approximation is n3 = 0.114213
An error of approximation of b by Ax for a random vector x in R^n is 40.4875
X =
   0.032693337120072   0.033395601743820
  -0.221874291154678  -0.136900271681156
   0.256229883897294   0.179917128116664

%E
A=magic(4);A=orth(A),b=rand(4,1)
A =
  -0.500000000000000   0.670820393249937   0.500000000000000
  -0.500000000000000  -0.223606797749979  -0.500000000000000
  -0.500000000000000   0.223606797749979  -0.500000000000000
  -0.500000000000000  -0.670820393249937   0.500000000000000
b =
   0.035711678574190
   0.849129305868777
   0.933993247757551
   0.678735154857773

X = solvemore(A,b)
The system is inconsistent - look for the least-squares solution
The solution of the normal equations is
  -1.248784693529146
  -0.412377106939306
  -0.534337860097182
The least-squares error of the approximation is n1 = 0.200713
A has orthonormal columns: an orthonormal basis for Col A is U=A
The projection of b onto Col A is
   0.035711678574190
   0.849129305868777
   0.634788147107582
   0.678735154857774
The least-squares solution by using the projection onto Col A is x4 =
  -1.099182143204161
  -0.479281401366106
  -0.384735309772198
The least-squares error of this approximation is n2 = 0.299205
The least-squares error of this approximation is n3 = 0.221896
An error of approximation of b by Ax for a random vector x in R^n is 2.50203
X =
  -1.248784693529146  -1.099182143204161
  -0.412377106939306  -0.479281401366106
  -0.534337860097182  -0.384735309772198
diary off










diary on
format compact
% Exercise 5
type polyplot

function []=polyplot(a,b,p)
x=(a:(b-a)/50:b)';
y=polyval(p,x);
plot(x,y)
end


type lstsqline

function c=lstsqline(x,y)
format rat,
x=x';
y=y';
a=x(1);
m=length(x);
b=x(m);
X=[x,ones(m,1)];
c=lscov(X,y);
c1=(inv(X'*X))*(X'*y)
c2=(X'*X)\(X'*y)
N=norm(y-X*c)
plot(x,y,'*'), hold
polyplot(a,b,c');
P=poly2sym(c)
end

x=[0 2 3 5 6]; y=[4 3 2 1 0]
y =
       4              3              2              1              0       
c=lstsqline(x,y)
c1 =
     -25/38    
      78/19    
c2 =
     -25/38    
      78/19    
N =
     514/1417  
Current plot held
P =
78/19 - (25*x)/38
c =
     -25/38    
      78/19    


 



















% Exercise 6

type lstsqpoly

function c=lstsqpoly(x,y,n)
format rat,
x=x';
y=y';
a=x(1);
m=length(x);
b=x(m);
% Create matrix X dependent on degree of polynomial
rows=length(x);
columns=n+1;
X=zeros(rows,columns);
for currentColumn=1:n+1
    X(:,currentColumn)=x.^(n+1-currentColumn);
end
c=lscov(X,y);
c1=(inv(X'*X))*(X'*y)
c2=(X'*X)\(X'*y)
N=norm(y-X*c)
plot(x,y,'*'), hold
polyplot(a,b,c');
P=poly2sym(c)
end

x =
       0              2              3              5              6       

y =
       4              3              2              1              0       

n=1
n =
       1       
c=lstsqpoly(x,y,n)
c1 =
     -25/38    
      78/19    
c2 =
     -25/38    
      78/19    
N =
     514/1417  
Current plot held
P =
78/19 - (25*x)/38
c =
     -25/38    
      78/19 
% c=c1=c2= [-25/38;78/19]
% Note that this result is equal to that of Exercise 5
 

n=2
n =
       2       
c=lstsqpoly(x,y,n)
c1 =
      -3/154   
     -83/154   
     309/77    
c2 =
      -3/154   
     -83/154   
     309/77    
N =
     703/2181  
Current plot held
P =
309/77 - (3*x^2)/154 - (83*x)/154
c =
      -3/154   
     -83/154   
     309/77
% c=c1=c2= [-3/154;-83/154;309/77]    

 

n=3
n =
       3       
c=lstsqpoly(x,y,n)
c1 =
      -1/228   
       5/266   
    -983/1596  
     535/133   
c2 =
      -1/228   
       5/266   
    -983/1596  
     535/133   
N =
     454/1425  
Current plot held
P =
- x^3/228 + (5417864213378195*x^2)/288230376151711744 - (983*x)/1596 + 535/133
c =
      -1/228   
       5/266   
    -983/1596  
     535/133 
% c=c1=c2= [-1/228;5/266;-983/1596;535/133]

   

n=4
n =
       4       
c=lstsqpoly(x,y,n)
c1 =
      -1/40    
      19/60    
     -51/40    
      59/60    
       4       
c2 =
      -1/40    
      19/60    
     -51/40    
      59/60    
       4       
N =
       1/82444248778916
Current plot held
P =
- x^4/40 + (19*x^3)/60 - (51*x^2)/40 + (59*x)/60 + 4
c =
      -1/40    
      19/60    
     -51/40    
      59/60
% c=c1=c2= [-1/40;19/60;-51/40;59/60]    
       4   
     

% Based on plots and norms, polynomial of degree 4 fits the data best. Visually, the plot corresponding to n=4 contains a regression line that passes through all points. The norm is 1/82444248778916 , very close to 0

diary off

