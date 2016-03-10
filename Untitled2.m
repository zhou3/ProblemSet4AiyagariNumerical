% PROGRAM NAME: ps4huggett.m
clear, clc
tic

% PARAMETERS
alpha=1/3;
beta = .99; %discount factor 
sigma = 2; % coefficient of risk aversion
delta=.25;
%rou=.5;
%sigma1=.2;
%ny=5;
%m=3;
%b = 0.5; % replacement ratio (unemployment benefits)
%y_s = [1, b]; % endowment in employment states
%PI = [.97 .03; .5 .5]; % transition matrix
lnZ=TAUCHEN1(5,.5,.2,3);
Z=exp(lnZ);
Z=Z';
PI=TAUCHEN2(5,.5,.2,3);
%PIstar=[ PI(2,1) / (1-PI(1,1) + PI(2,1)),...
    %1 - ( PI(2,1) / (1-PI(1,1) + PI(2,1)) ) ];
 %invariant distribution for transition matrix 
invPI=invdist(PI);
% aggregate labor supply 
N=Z*invPI';
% ASSET VECTOR(endogenous variable)
a_lo = 0; %lower bound of grid points
a_hi = 6;% guess upper bound of grid points
%num_a = 701;
num_a = 500;
% Guess the K
K=1.545;
%compute factor price given our guess K:
r=alpha*(N^(1-alpha))*(K^(alpha-1))+(1-delta);
w=(K^alpha)*(1-alpha)*(N^(-alpha));
a = linspace(a_lo, a_hi, num_a); % asset (row) vector


% INITIAL GUESS FOR q
%q_min = 0.98;
%q_max = 1;
%q_guess = (q_min + q_max) / 2;

% ITERATE OVER ASSET PRICES
%aggsav = 1 ;
%while abs(aggsav) >= 0.01 ;

% CURRENT RETURN (UTILITY) FUNCTION
cons = bsxfun(@minus, r*a',  a);
cons = bsxfun(@plus, cons, permute(Z*w, [1 3 2]));
ret = (cons .^ (1-sigma)) ./ (1 - sigma); % current period utility
ret(cons<0) = -Inf;

% INITIAL VALUE FUNCTION GUESS
v_guess = zeros(5, num_a);

% VALUE FUNCTION ITERATION
v_tol = 1;
while v_tol >.0001;
   % CONSTRUCT TOTAL RETURN FUNCTION
   v_mat = ret + beta * ...
       repmat(permute(PI * v_guess, [3 2 1]), [num_a 1 1]);
   
   % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
   [vfn, pol_indx] = max(v_mat, [], 2);
   vfn = permute(vfn, [3 1 2]);
   
   v_tol = abs(max(v_guess(:) - vfn(:)));
   
   v_guess = vfn; %update value functions
end;

% KEEP DECSISION RULE
pol_indx = permute(pol_indx, [3 1 2]);
pol_fn = a(pol_indx);

% SET UP INITITAL DISTRIBUTION
%Mu = zeros(1,num_a);

%Mu(1, 4) = 1; % initial guess: everyone employed, 0 assets
Mu_guess = ones(5, num_a); 
%alternative initial guess: same mass in all states
 Mu = Mu_guess / sum(Mu_guess(:)); % normalize total mass to 1

% ITERATE OVER DISTRIBUTIONS
% way 1: loop over all non-zeros states
mu_tol = 1;
while mu_tol > 1e-08
    [emp_ind, a_ind] = find(Mu > 0); % find non-zero indices
    
    MuNew = zeros(size(Mu));
    for ii = 1:length(emp_ind)
        apr_ind = pol_indx(emp_ind(ii), a_ind(ii)); 
        MuNew(:, apr_ind) = MuNew(:, apr_ind) + ...
            (PI(emp_ind(ii), :) * Mu(emp_ind(ii), a_ind(ii)) )';
    end

    mu_tol = max(abs(MuNew(:) - Mu(:)));
    
    Mu = MuNew ;
end

% % way 2: use transition matrices
% T_tilde = zeros(num_a, num_a, 2, 2);
% % set up matrices:
% for from_a = 1:num_a
%     for from_s = 1:2
%         T_tilde(from_a, pol_indx(from_s, from_a), from_s, :) = PI(from_s,:);
%     end
% end
% % transition:
% mu_tol = 1;
% while mu_tol > 1e-08
%     MuNew = zeros(size(Mu));
%     for from_s = 1:2
%         for to_s = 1:2
%             MuNew(to_s,:) = MuNew(to_s,:) + Mu(from_s,:) * T_tilde(:,:,from_s,to_s);
%         end
%     end
%     mu_tol = max(abs(Mu(:) - MuNew(:)));
%     Mu = MuNew;
% end


% CHECK AGGREGATE DEMAND
aggsav = sum( pol_fn(:) .* Mu(:) ); % Aggregate future assets

%if aggsav > 0 ;
 %   q_min = q_guess ;
%end ;
%if aggsav < 0;
 %   q_max = q_guess ;
%end ;

%display (['q = ', num2str(q_guess)])
%display (['Aggregate desired wealth = ', num2str(aggsav)]);
%display (['New qmin is ', num2str(q_min), ', new qmax is ', num2str(q_max)]);
%display (['New q is ', num2str((q_max + q_min)/2)]);

%q_guess = (q_max + q_min)/2 ;

%display (' ') ;

%end
%after several calibration, we know K should equail 1.54.
%Mu is the steady state distribution
%in SS ,r=1.005, r(cm)=1/beta=1.0101;
plot(a,pol_fn(1,:));
plot(a,pol_fn(2,:));
plot(a,pol_fn(3,:));
plot(a,pol_fn(4,:));
plot(a,pol_fn(5,:));

%agg_wealth = sum(Mu,1) * a' + y_s * sum(Mu,2); 
% wealth is asset holdings plus incomes
agg_wealth = sum(Mu,1) * a';
%wealth_dist = [[Mu(1,:), Mu(2,:)]; [a + y_s(1), a + y_s(2)]]';
wealth_dist = [[Mu(1,:), Mu(2,:),Mu(3,:), Mu(4,:),Mu(5,:)]; a,a,a,a,a]';
[~, ordr] = sort(wealth_dist(:,2), 1);
wealth_dist = wealth_dist(ordr,:);
% compute gini
pct_dist = cumsum( (wealth_dist(:,2) ./ agg_wealth) .* wealth_dist(:,1) );
gini = 1 - sum( ([0; pct_dist(1:end-1)] + pct_dist) .* wealth_dist(:,1) );
x= linspace(0, 1, 2500);
cumsum( (wealth_dist(:,1)).*(wealth_dist(:,2))/ (agg_wealth));
plot(x,pct_dist);
hold on 
plot(x,x');
hold off