% PARAMETERS
% use linear interpolation
alpha=1/3;
beta = .99; %discount factor 
sigma = 2; % coefficient of risk aversion
delta=.025;
% original method to compute the value function:
% get the state variable Z 
lnZ=TAUCHEN1(5,.5,.2,3);
Z=exp(lnZ);
PI=TAUCHEN2(5,.5,.2,3);
invPI=invdist(PI);
% aggregate labor supply 
N=invPI*Z;
% get transition matrix
%PI=TAUCHEN2(5,.5,.2,3);
a_lo = 0; %lower bound of grid points
a_hi = 40;% guess upper bound of grid points
%num_a = 701;
num_a = 50;
a=linspace(a_lo,a_hi,num_a);
%initialize the value function
%vfn=zeros(500,500,5);
%ret_fn=(cons .^ (1-sigma)) ./ (1 - sigma);
% compute interest rate
K=29;
%as we already solve 
r=alpha*(N^(1-alpha))*(K^(alpha-1))+(1-delta);
%r=alpha*(N^(1-alpha))*(K^(alpha-1))+(1-delta);
w=(K^alpha)*(1-alpha)*(N^(-alpha));
%ret_fn=(1/(1 - sigma))*(cons .^ (1-sigma));
% compute and create the matrix cons:
%cons=zeros(500,500,5);
%repmat(Z*w,[1,num_a,num_a]);
ZW=permute(repmat(Z*w,[1,num_a,num_a]),[2 3 1]);
%for n=1:5
%WZ=zeros(num_a,num_a,5);
%WZ(:,:,n)=Z(n,1)*w;
%end 
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
v_tol = 1;
v_guess=zeros(num_a,5);
tic 
while v_tol>.0001;
%vfn=zeros(500,5);
%vfn_mat=zeros(500,500,5);

%get the expected vfn
exp_vfn=v_guess*PI';
exp_vfn=permute(exp_vfn,[1 3 2]);
exp_finvfn=repmat(exp_vfn,[1 num_a 1]);
exp_finvfn=permute(exp_finvfn,[2 1 3]);
%for ii=1:10
 %   for jj=1:5
  %  exp_finvfn(:,ii,jj)=exp_vfn(ii,jj);
   % end 
%end 

%exp_finvfn=repmat(exp_vfn,[]);
v_mat=ret_fn+beta*exp_finvfn;
 [vfn, pol_indx] = max(v_mat, [], 2);
%haha=max(hehe,[],2);
vfn=permute(vfn,[1 3 2]);
% interpolate the original value function
v_tol=abs(max(v_guess(:)-vfn(:)));
%abs(max(v_guess(:) - vfn(:)));
v_guess=vfn;
end

% interpolate the original value function
num_a2=500;
aa=linspace(a_lo,a_hi,num_a2);
%define query point
v_guess=permute(v_guess,[1 3 2]);
v_betterguess=interp1(a',v_guess,aa');
v_betterguess=permute(v_betterguess,[1 3 2]);
vv_tol=1;
% New cons
ZW1=permute(repmat(Z*w,[1,num_a2,num_a2]),[2 3 1]);
RA1=zeros(num_a2,num_a2,5);
for i=1:num_a2
RA1(i,:,:)=aa(i)*r;
end 
A_pri1=zeros(num_a2,num_a2,5);
for j=1:num_a2
A_pri1(:,j,:)=aa(j);
end 
cons=ZW1+RA1-A_pri1;
%ret_fn=(1/(1 - sigma))*(cons .^ (1-sigma));
ret_fn1 = (cons .^ (1-sigma)) ./ (1 - sigma);
ret_fn1(cons<0) = -Inf;

while vv_tol>.0001;
%vfn=zeros(500,5);
%vfn_mat=zeros(500,500,5);

%get the expected vfn
exp_vfn=v_betterguess*PI';
exp_vfn=permute(exp_vfn,[1 3 2]);
exp_finvfn=repmat(exp_vfn,[1 num_a2 1]);
exp_finvfn=permute(exp_finvfn,[2 1 3]);
%for ii=1:10
 %   for jj=1:5
  %  exp_finvfn(:,ii,jj)=exp_vfn(ii,jj);
   % end 
%end 

%exp_finvfn=repmat(exp_vfn,[]);
v_mat=ret_fn1+beta*exp_finvfn;
 [vfn1, pol_indx1] = max(v_mat, [], 2);
%haha=max(hehe,[],2);
vfn1=permute(vfn1,[1 3 2]);
% interpolate the original value function
vv_tol=abs(max(v_betterguess(:)-vfn1(:)));
%abs(max(v_guess(:) - vfn(:)));
v_betterguess=vfn1;
end
toc
%x = 0:pi/4:2*pi;
%v = sin(x);
%plot(x,v,'.')
%xq = 0:pi/16:2*pi;
% check whether distribution stay same or not
pol_indx1 = permute(pol_indx1, [3 1 2]);
%pol_fn1 = aa'(pol_indx1);

Mu_guess = ones(5, num_a2); 
%alternative initial guess: same mass in all states
 Mu = Mu_guess / sum(Mu_guess(:)); % normalize total mass to 1

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1;
while mu_tol > 1e-08
    [emp_ind, a_ind] = find(Mu > 0); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx1(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
    end

    mu_tol = max(abs(MuNew(:) - Mu(:)));
    
    Mu = MuNew ;
end
%Mu=permute(Mu,[2 3 1]);
plot(aa',Mu);
