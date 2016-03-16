alpha=1/3;
beta = .99; %discount factor 
sigma = 2; % coefficient of risk aversion
delta=.025;
lnZ=TAUCHEN1(5,.5,.2,3);
Z=exp(lnZ);
PI=TAUCHEN2(5,.5,.2,3);
invPI=invdist(PI);
N=invPI*Z;
% get transition matrix
%PI=TAUCHEN2(5,.5,.2,3);
a_lo = 0; %lower bound of grid points
a_hi = 40;% guess upper bound of grid points
%num_a = 701;
num_a = 500;
a=linspace(a_lo,a_hi,num_a);
K=29;
r=alpha*(N^(1-alpha))*(K^(alpha-1))+(1-delta);
w=(K^alpha)*(1-alpha)*(N^(-alpha));
ZW=permute(repmat(Z*w,[1,num_a,num_a]),[2 3 1]);
RA=zeros(num_a,num_a,5);
for i=1:num_a
RA(i,:,:)=a(i)*r;
end 
A_pri=zeros(num_a,num_a,5);
for j=1:num_a
A_pri(:,j,:)=a(j);
end 
cons=ZW+RA-A_pri;
%ret_fn=(1/(1 - sigma))*(cons .^ (1-sigma));
ret_fn = (cons .^ (1-sigma)) ./ (1 - sigma);
ret_fn(cons<0) = -Inf;
%initialize value_fun:
v_tol = 100;
v_guess=zeros(num_a,5);
%tic 
while v_tol>.0001;
exp_vfn=v_guess*PI';
exp_vfn=permute(exp_vfn,[1 3 2]);
exp_finvfn=repmat(exp_vfn,[1 num_a 1]);
exp_finvfn=permute(exp_finvfn,[2 1 3]);

%w=vfn;
%for j=1:30;
v_mat=ret_fn+beta*exp_finvfn;
 [vfn, pol_indx] = max(v_mat, [], 2);
%haha=max(hehe,[],2);
vfn=permute(vfn,[1 3 2]);% get 500*5
vfn=vfn';% get 5*500;
vfn=reshape(vfn,[5*num_a 1]);
newcons=zeros(num_a,5);
for ii=1:num_a
    for jj=1:5
newcons(ii,jj)=cons(ii,pol_indx(ii,1,jj),jj);
    end 
end 
newcons=newcons';
retnew = (newcons .^ (1-sigma)) ./ (1 - sigma);
retnew(newcons<0) = -Inf;
retnew=reshape(retnew,[5*num_a 1]);
% define matrix Q
pol_indx1=permute(pol_indx,[3 1 2]);
QQ=makeQmatrix(pol_indx1,PI);
v_new=zeros(5*num_a,1,31);
for k=1:30;
    v_new(:,:,1)=retnew+beta*QQ*vfn;
    v_new(:,:,k+1)=retnew+beta*QQ*v_new(:,:,k);
end 
v_fin=v_new(:,:,31);
% interpolate the original value function
v_guess=v_guess';
v_tol=abs(max(v_guess(:)-v_fin(:)));
%abs(max(v_guess(:) - vfn(:)));
v_guess=reshape(v_fin,[5 num_a]);
v_guess=v_guess';
end
    