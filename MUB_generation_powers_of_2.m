clc
clear all
tic
n=12;
N=2^n; % MUB for C^N 
Nb=64; % Nb<=N (excluding eye) .........number of bases of C^N to be calculated 
Nb=min(N,Nb);
msbfirst = false;
bin=int2bit(0:2^n-1,n,msbfirst)';
binG=bin*bin';
sign= 1-2*(mod(binG,2)==1);
%% list of primitive polynomials of degree degree n
h{1}=[1 1];
h{2}=[1 1 1];
h{3}=[1 1 0 1];
h{4}=[1 1 0 0 1];
h{5}=[1 0 1 0 0 1];
h{6}=[1 1 0 0 0 0 1];
h{7}=[1 1 0 0 0 0 0 1];
h{8}=[1 0 1 1 1 0 0 0 1];
h{9}=[1 0 0 0 1 0 0 0 0 1];
h{10}=[1 0 0 1 0 0 0 0 0 0 1];
h{11}=[1 0 1 0 0 0 0 0 0 0 0 1];
h{12}=[1 1 0 0 1 0 1 0 0 0 0 0 1];
%% calculation for alpha matrices
alpha=zeros(n,n,n);
a=h{n};
for i=1:n
    for j=1:n
        b=zeros(1,i+j-1);
        b(i+j-1)=1;
        [q,r]=gfdeconv(b,a);
        c=find(r);
        for m=1:length(c)
            alpha(i,j,c(m))=1;
        end
    end
end
%% calculation of bases
b=zeros(N,N);
B=[];
for r=0:Nb-1
    r_binary=bin(r+1,:);
    r_alpha=zeros(n);
    for i=1:n
        r_alpha=r_alpha+r_binary(i)*alpha(:,:,i);
    end
    for l=0:N-1
        b(:,l+1)=(1i^(bin(l+1,:)*r_alpha*bin(l+1,:)'));
    end
    B=[B (b.*sign).'];
end
B=B/sqrt(N);
toc
% G=unique(abs(B'*B))
string="first "+string(r+1)+" MUBs of dimension 2^"+string(n);
txt="MUB_2^"+string(n)+".mat";
save('txt','B','string')
