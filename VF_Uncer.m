clear
close all
clc
%%%% Set up parameters
alpha = 0.35;
beta = 0.99;
delta = 0.025;
sigma = 2;

%%%% Set up discretized state space
k_min = 0;
k_max = 45;
num_k = 1000; % number of points in the grid for k

k = linspace(k_min, k_max, num_k);

k_mat = repmat(k', [1 num_k]); % this will be useful in a bit
%Create A_t
A1=zeros(999,1);
A2=zeros(999,1);
A1(1)=1.1;
A2(1)=0.678;
for i=1:size(A1,1)
    if A1(i)>1
         if binornd(1,0.977)>0.5
             A1(i+1)=1.1;
         else
             A1(i+1)=0.678;
         end
    else
        if binornd(1,0.926) >0.5
            A1(i+1) =0.678;
        else
            A1(i+1) =1.1;
        end
    end
end
%for i=1:size(A2,1)
    %if A2(i)<1
       if binornd(1,0.926)>0.5
           A2(i+1)=0.678;
       else
           A1(i+1)=1.1;
       end
    else
        if binornd(1,0.977)>0.5
            A2(i+1)=1.1;
        else
            A2(i+1)=0.678;
        end
    end
end

%%%% Set up consumption and return function
% 1st dim(rows): k today, 2nd dim (cols): k' chosen for tomorrow
cons1= A1.*k_mat .^ alpha + (1 - delta) * k_mat - k_mat'; 
cons2= A2.*k_mat .^ alpha + (1 - delta) * k_mat - k_mat';

ret1 = cons1 .^ (1 - sigma) / (1 - sigma);
ret2 = cons2 .^ (1 - sigma) / (1 - sigma);% return function
% negative consumption is not possible -> make it irrelevant by assigning
% it very large negative utility
ret1(cons1 < 0) = -Inf;
ret2(cons2 < 0) = -Inf;

%%%% Iteration
dis = 1; dis2=1; tol = 1e-06; % tolerance for stopping 
v_guess1 = zeros(1, num_k);
v_guess2 = zeros(1, num_k);
while dis > tol  
    % compute the utility value for all possible combinations of k and k':
    value_mat1 = ret1 + beta * repmat(v_guess1, [num_k 1]);
    
    % find the optimal k' for every k:
    [vfn1, pol_indx1] = max(value_mat1, [], 2);
   
    vfn1 = vfn1';
    
    
    % what is the distance between current guess and value function
    dis = max(abs(vfn1 - v_guess1));
  
    
    % if distance is larger than tolerance, update current guess and
    % continue, otherwise exit the loop
    v_guess1 = vfn1;
 
end
while dis2 > tol
    value_mat2 = ret2 + beta * repmat(v_guess2, [num_k 1]);
     [vfn2, pol_indx2] = max(value_mat2, [], 2);vfn2 = vfn2';
       dis2= max(abs(vfn2 - v_guess2));
        v_gusee2 = vfn2;
end
     
g1 = k(pol_indx1); % policy function
g2 = k(pol_indx2);

plot(k,vfn1)
figure
plot(k,vfn2)
figure
plot(k,g1)
figure
plot(k,g2)
