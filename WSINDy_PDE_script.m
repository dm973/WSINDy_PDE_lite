%% load data

load('burgers.mat')

coarsen_data = [[0 1 1];[0 1 1]];
Shift = [];
sigma_NR = 0.2;
rng('shuffle');
noise_dist = 0;
noise_alg = 0;
rng_seed = rng().Seed; 
rng(rng_seed);
toggle_disp = 1;

[U_obs,xs_obs] = coarsen_data_fcn(U_exact,xs,coarsen_data,Shift);
[U_obs,noise,snr,sigma] = gen_noise(U_obs,sigma_NR,noise_dist,noise_alg,toggle_disp);

%% specify library

max_dx = 5;
max_dt = 1;
polys = [0:5];
trigs = [];
use_all_dt = 0;
use_cross_dx = 0;
custom_add = [];
custom_remove = {};
toggle_comb = 0;

nonautargs = {}; 
convargs = {};
drifttags = {};

n = length(U_obs); 
dims = size(U_obs{1});
dim = length(dims);
[tags_pde,lib_list,~,lhs_ind,max_dx,max_dt,polys,customf,customconv]  = get_lib_tags(n,dim,lhs,max_dx,max_dt,polys,trigs,use_cross_dx,use_all_dt,custom_remove,custom_add,drifttags,convargs,true_nz_weights);

%% Set Discretization

m_x = 15;
p_x = 9; 
phifun_x = @(x) (1-x.^2).^p_x;

m_t = 15;
p_t = 9;
phifun_t = @(x) (1-x.^2).^p_t;

sm_x = 5;
sm_t = 5;

[sub_inds,ss] = get_subinds(xs_obs,m_x,m_t,sm_x,sm_t);

Cfs_x = phi_weights(phifun_x,m_x,max_dx);
Cfs_t = phi_weights(phifun_t,m_t,max_dt);

%% Set scales

scales = 2;

ms = [m_x m_t]; ps = [p_x p_t]; phi_class = {1,1};
[scales,M_full] = get_scales(U_obs,scales,polys,ps,ms,max_dx,max_dt,lib_list,xs_obs,customf,customconv,phi_class);

%% Build Linear system

[Theta_pdx,libtree] = get_lib_columns(n,lib_list,U_obs,Cfs_x,Cfs_t,m_x,m_t,sub_inds,dim,scales,xs_obs,customf,customconv);
[G,b,M] = get_linear_system(Theta_pdx,lhs_ind,M_full);

%% Solve sparse regression problem

lambda = 10.^(linspace(-4, 0, 50));
gamma = 0;
alpha = 0.5;
maxits = 20;

[W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq(lambda,gamma,G,b,M,maxits,alpha);

%% display learned PDE and algorithm info

tags_pde_G = tags_pde(~ismember(1:length(tags_pde),lhs_ind));
lib_list_G = lib_list(~ismember(1:size(lib_list,1),lhs_ind),:);
axi = tags2axi(true_nz_weights,lib_list);
axi = axi(~ismember(1:length(tags_pde),lhs_ind),:);

lambda_hat_ind=find(lossvals(min(end,4),:)>0,1,'last');
lambda_hat = lossvals(min(end,4),lambda_hat_ind);

print_loc=1;
[m,n] = size(W);
str_wsindy = cell(n,1);
for k=1:n
    str_wsindy{k} = print_pde(W(:,k),tags_pde_G,tags_pde{lhs_ind(k)});
end
if ~isempty(axi)
    Tps = tpscore({W},axi);
else
    Tps=NaN;
end
dW = wnorm({W},axi,Inf);

if ~isequal(print_loc,0)
    if ~isequal(print_loc,1)
        print_loc = fopen(print_loc,'a');
    end

    for k=1:n
        fprintf(print_loc,['\nRecovered PDE: ',str_wsindy{k}]);
        fprintf(print_loc,'\nRelative Res: ||b-G*W||_2/||b||_2 = %.2e',norm(resid(:,k)));
        if ~isempty(dW)
            fprintf(print_loc,'\nMax Weight Error: max|W-W_{true}| = %.2e\n', dW(k));
        end
    end
    fprintf(print_loc,'TP Score = %1.2f\n', Tps);
    fprintf(print_loc,'      \n');
    fprintf(print_loc,'polys = ');
    fprintf(print_loc,'%u ',polys);
    fprintf(print_loc,'\ntrigs = ');
    fprintf(print_loc,'%u ',trigs);
    fprintf(print_loc,'\nMax derivs [t x] = ');
    fprintf(print_loc,'%u ',[max_dt max_dx]);
    fprintf(print_loc,'\n[m_x m_t] = ');
    fprintf(print_loc,'%u ',[m_x m_t]);
    fprintf(print_loc,'\n[s_x s_t] = ');
    fprintf(print_loc,'%u ',[diff(sub_inds{1}(1:min(length(sub_inds{1}),2))') diff(sub_inds{end}(1:min(length(sub_inds{end}),2))')]);
    fprintf(print_loc,'\n[p_x p_t] = ');
    pps=zeros(1,2);
    ps=[p_x p_t];
    for i=1:2
        if isequal(phi_class{i},1)
            pps(i)=ps(i);
        elseif isequal(phi_class{i},2)
            pps(i)=1/ps(i);
        else
            pps(i)=NaN;
        end
    end

    fprintf(print_loc,'%u ',pps);
    fprintf(print_loc,'\n scales = ');
    fprintf(print_loc,'%.2e ',scales) ;
    fprintf(print_loc,'\n      \n');
    fprintf(print_loc,'Size of dataset = ');
    fprintf(print_loc,'%u ',dims);
    fprintf(print_loc,'\nSize G = ');
    fprintf(print_loc,'%u ',size(G));
    if gamma >0
        fprintf(print_loc,'\nCond G = %.2e',cond([G;gamma*norm(G)*eye(m)]));
    else
        fprintf(print_loc,'\nCond G = %.2e',cond(G));
    end
    fprintf(print_loc,'\n[lambda_hat gamma] = ');
    fprintf(print_loc,'%.3e ',[lambda_hat gamma*norm(G)]);
    fprintf(print_loc,'\n[sigma_NR sigma] = ');
    fprintf(print_loc,'%.3e ',[sigma_NR sigma]);
    fprintf(print_loc,'\n      \n');
    fprintf(print_loc,'STLS its = ');
    fprintf(print_loc,'%u ',its_all);
    fprintf(print_loc,'\n ');

    if ~all(print_loc==1)
        fclose(print_loc);
    end
end

%% plot solution

for nn=1:n
    figure(nn)
    surf(xs_obs{1},xs_obs{2},U_obs{nn}','edgecolor','none')
    view([0 90])
    legend({['u',num2str(nn)]})
    colorbar
    xlabel('x')
    ylabel('t')
end

function [U_obs,noise,snr,sigma] = gen_noise(U_exact,sigma_NR,noise_dist,noise_alg,toggle_disp)

    n = length(U_exact);

    if sigma_NR~=0
        U_obs = cell(n,1);
        stdvs = zeros(1,n);
        noise = cell(1,n);
        snr = zeros(1,n);
        for k=1:n
            stdvs(k) = rms(U_exact{k}(:))^2;
        end
        for j=1:n
            [U_obs{j},noise{j},snr(j),sigma] = add_noise(U_exact{j},stdvs(j),sigma_NR,noise_dist,noise_alg);
            if toggle_disp
                disp(['[SNR sigma] = ',num2str([snr(j) sigma])]); 
            end
        end
    else
        U_obs = U_exact;
        noise = [];
        snr = 0;
        sigma = 0;
    end

end


function [U_obs,xs_obs] = coarsen_data_fcn(U_exact,xs,coarsen_data,Shift)

    n = length(U_exact);
    dims = size(U_exact{1});
    dim = length(dims);
    
    U_obs = U_exact;
    xs_obs = xs;
    
    if ~isempty(coarsen_data)
        inds = cell(1,dim);
        for j=1:dim
            inds{j} = 1+floor(coarsen_data(j,1)*dims(j)):ceil(coarsen_data(j,2)):ceil(coarsen_data(j,3)*dims(j));
            xs_obs{j} = xs{j}(inds{j});
        end
        for j=1:n
            U_obs{j} = U_exact{j}(inds{:});
        end
    end
    
    for j=1:min(length(Shift),dim)
        xs_obs{j} = xs_obs{j}+Shift(j);
    end

end

function [sub_inds,ss] = get_subinds(xs_obs,m_x,m_t,sm_x,sm_t)
    dims = cellfun(@(x)length(x),xs_obs);
    dim = length(xs_obs);
    s_x = max(floor(m_x/sm_x),1);max(floor(length(xs_obs{1})/25),1);
    s_t = max(floor(m_t/sm_t),1);max(floor(length(xs_obs{end})/25),1);

    sub_inds = cell(1,dim);
    ss = [repmat(s_x,1,dim-1) s_t];
    mm = [repmat(m_x,1,dim-1) m_t];
    for j=1:dim
        N = dims(j);
        m = mm(j);
        s = ss(j);
        end_pt=N-2*m;
        sub_inds{j} = 1:s:end_pt;
    end
end

function [W,resid,its_all,lossvals,thrs_EL] = wsindy_pde_RGLS_seq(lambdas,gamma,G,b,M,maxits,alpha)
    
    [~,m] = size(G);
    [~,num_eq] = size(b);
    
    W_ls = [G;gamma*eye(m)] \ [b;zeros(m,num_eq)];
    GW_ls = norm(G*W_ls);
    
    proj_cost = [];
    overfit_cost = [];
    lossvals = [];
    
    if isempty(lambdas)
        lam_max = max(max(abs(G'*b),[],2)./vecnorm(G).^2');
        lam_min = min(vecnorm(G*W_ls))/size(G,2)/max(vecnorm(G));
        lambdas = 10.^linspace(log10(lam_min), log10(lam_max),50);
    end
    
    if and(length(lambdas)==1,all(lambdas<0))
        lam_max = max(max(abs(G'*b),[],2)./vecnorm(G).^2');
        lam_min = min(vecnorm(G*W_ls))/size(G,2)/max(vecnorm(G));
        lambdas = 10.^linspace(log10(lam_min), log10(lam_max),-lambdas);
    end
    
    W = zeros(m,num_eq);
    for l=1:length(lambdas)
        lambda = lambdas(l);
        for k=1:num_eq
            if isempty(M)
                [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, [], maxits);
            else
                [W(:,k),~,~] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,k), maxits);
                W(:,k) = W(:,k)./M(:,k);
            end
        end
        proj_cost = [proj_cost 2*alpha*norm(G*(W-W_ls))/GW_ls];
        overfit_cost = [overfit_cost 2*(1-alpha)*length(find(W~=0))/length(find(W_ls~=0))];
        lossvals = [lossvals proj_cost(end) + overfit_cost(end)];
    end
    
    l = find(lossvals == min(lossvals),1);
    
    lambda = lambdas(l);
    its_all = zeros(num_eq,1);
    
    resid = b*0;
    for k=1:num_eq
        if isempty(M)
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, []);
            resid(:,k) = (b(:,k) - G*W(:,k))/norm(b(:,k)); 
        else
            [W(:,k),its,thrs_EL] = sparsifyDynamics(G, b(:,k), lambda, gamma, M(:,k), maxits);
            resid(:,k) = (b(:,k) - G*(W(:,k)./M(:,k)))/norm(b(:,k)); 
        end
        its_all(k) = its;
    end
    lossvals = [lossvals;lambdas; [[lossvals(1:l);lambdas(1:l)] zeros(2,length(lambdas)-l)]; proj_cost; overfit_cost];
end


function [Xi,its,thrs_EL] = sparsifyDynamics(Theta,dXdt,lambda,gamma,M,maxits)
% Copyright 2015, All Rights Reserved
% Code by Steven L. Brunton
% For Paper, "Discovering Governing Equations from Data: 
%        Sparse Identification of Nonlinear Dynamical Systems"
% by S. L. Brunton, J. L. Proctor, and J. N. Kutz
%
% modified by Daniel A. Messenger, 2020 to prevent return of zero vector
% and include regularization
%
% compute Sparse regression: sequential least squares

n = min(size(dXdt));
nn = size(Theta,2);

if ~exist('maxits','var')
    maxits= nn;
end

if  gamma ~= 0
    Theta_reg = [Theta;gamma*norm(Theta)*eye(nn)];
    dXdt_reg = [dXdt;zeros(nn,n)];
else
    Theta_reg = Theta;
    dXdt_reg = dXdt;
end

Xi = Theta_reg \ dXdt_reg;  % initial guess: Least-squares
if ~isempty(M)
    Xi = M.*Xi;
end

if isempty(M)
    thrs_EL = [];
else
    bnds = norm(dXdt_reg)./vecnorm(Theta_reg)'.*M; 
    LBs = lambda*max(1,bnds);
    UBs = 1/lambda*min(1,bnds);
    thrs_EL = [LBs bnds UBs];
end

smallinds = 0*Xi;
its = 0;
while its < maxits
    if ~isempty(M)
        smallinds_new = or(abs(Xi)<LBs,abs(Xi)>UBs);
        if or(length(find(smallinds_new))==nn,all(smallinds_new(:)==smallinds(:)))
            return
        else
            smallinds = smallinds_new;
            if  gamma ~= 0
                Theta_reg = [Theta;gamma*norm(Theta(:,~smallinds))*eye(nn)];
            end
            Xi(smallinds)=0;    
            for ind=1:n
                Xi(~smallinds,ind) = M(~smallinds).*(Theta_reg(:,~smallinds)\dXdt_reg(:,ind));
            end
        end
    else
        smallinds_new = (abs(Xi)<lambda);
        if or(length(find(smallinds_new))==nn,all(smallinds_new(:)==smallinds(:)))
            return
        else
            smallinds = smallinds_new;
            if  gamma ~= 0
                Theta_reg = [Theta;gamma*norm(Theta(:,~smallinds))*eye(nn)];
            end
            Xi(smallinds)=0;
            for ind = 1:n        
                biginds = ~smallinds(:,ind);
                Xi(biginds,ind) = Theta_reg(:,biginds)\dXdt_reg(:,ind);
            end
        end
    end
    its = its + 1;
end
end

function Cfs = phi_weights(phifun,m,maxd)
    x = linspace(-1,1,2*m+1);
    Cfs = zeros(maxd+1,2*m+1);
    syms y;
    f = @(y)phifun(y);
    for j=1:maxd+1
        Df = matlabFunction(diff(f(y),j-1));
        Cfs(j,:) = fillmissing(Df(x),'constant',Df(eps));
        inds = find(isinf(abs(Cfs(j,:))));
        for k=1:length(inds)
            Cfs(j,inds(k)) = Df(x(inds(k))+eps);
        end
    end
    Cfs = Cfs/norm(Cfs(1,:),1);
end

function [G,b,M] = get_linear_system(Theta_pdx,lhs_ind,M_full)

    num_eq = length(lhs_ind);
    [K,m] = size(Theta_pdx);
    G = Theta_pdx(:,~ismember(1:m,lhs_ind));
    b = zeros(K,num_eq);
    M = [];
    for k=1:num_eq
        b(:,k) = Theta_pdx(:,lhs_ind(k));
        if ~isempty(M_full)
            M = [M M_full(~ismember(1:m,lhs_ind))/M_full(lhs_ind(k))];
        end
    end
end
