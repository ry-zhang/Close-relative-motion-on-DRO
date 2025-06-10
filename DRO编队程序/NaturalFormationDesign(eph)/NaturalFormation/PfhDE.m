% 差分进化算法，该程序适应于求解的优化问题类型为：百个变量以内、连续函数和无约束
% v0：朱小龙(2020/07/07)
%     创建并撰写核心程序
% v1：杨驰航(2022/02/18),
%     将单独的max_nfes设置扩展至opts_DE的五个可调节属性；
%     去掉输入参数problem_size，改为以lb的size计算

%% PfhDE算法优化主函数

function [bsf_fit, bsf_sol] = PfhDE(fx,lb,ub,opts_DE)
%% 输入输出说明
% fx:               需要优化的目标函数
% problem_size:     优化变量个数
% lb:               优化变量最小值,1*problem_size（优化变量个数）数组
% ub:               优化变量最大值,1*problem_size数组
% bsf_fit:          每一代的最优值，G_Max（最大迭代代数）*1数组
% bsf_sol:          每一代的最优解，G_Max*problem_size矩阵
% opts_DE.Max_nfes: maximum fun-count (integer, default as 10000*problem_size)
%        .Display:  display the optimization process or not (0/1)
%        .Resize:   resize the population or not when nfes>max_nefs/2 (0/1)
%        .UseParallel: apply parallel computing or not (0/1)
%        .PopSize:  initial population size (integer)
% 
%% 实例程序
%  fx = @(x)sum(x.^2);
%  problem_size = 2;
%  lb = -10*ones(1,problem_size);
%  ub = 10*ones(1,problem_size);
%  opts_DE.Max_nfes = problem_size*10000;
%  opts_DE.Display = 1; opts_DE.Resize = 0;
%  opts_DE.UseParallel = 1; 
%  [bsf_fit, bsf_sol] = PfhDE(fx,lb,ub,opts_DE);

%% opts参数传递
max_nfes = opts_DE.Max_nfes; % maximum fun-count
isdisplay = opts_DE.Display; % display the process or not 
isresize = opts_DE.Resize; % resize the popuolation or not when nfes>max_nefs/2
UseParallel = opts_DE.UseParallel;

%% DE本身参数
problem_size = max(size(lb));
freq_inti = 0.5;
pb = 0.4;
ps = .5;
S.Ndim = problem_size;
S.Lband = lb;
S.Uband = ub;
lu = [lb;ub];
GenMaxSelected = 250; %% For local search
%% DE本身参数

%% parameter settings for L-SHADE
G_Max = G_Max_calc(problem_size,max_nfes);
p_best_rate = 0.11;    %0.11
arc_rate = 1.4;
memory_size = 5;

if isfield(opts_DE,'PopSize')
    pop_size = opts_DE.PopSize;
else
    pop_size = 18 * problem_size;   %18*D
end

SEL = round(ps*pop_size);

max_pop_size = pop_size;
min_pop_size = 4.0;

run_funcvals = [];

nfes = 0;
%% Initialize the main population
popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .* ...
    (repmat(lu(2, :) - lu(1, :), pop_size, 1));
pop = popold; % the old population becomes the current population

fitness = zeros(1,pop_size);
if UseParallel == 1
    parfor m = 1:pop_size
        fitness(m) = fx(pop(m,:));
    end
else
    for m = 1:pop_size
        fitness(m) = fx(pop(m,:));
    end
end
fitness = fitness';

%%% Initialize LS population
flag_LS = false;
counter = 0;
popsize_LS = 10;

%%% Initialize LS population to re-start them
popLS = repmat(lu(1, :), popsize_LS, 1) + rand(popsize_LS, problem_size) .* ...
    (repmat(lu(2, :) - lu(1, :), popsize_LS, 1));
fitness_LS = zeros(1,popsize_LS);
for m = 1:popsize_LS
    
    fitness_LS(m) = fx(popLS(m,:));
    
end

fitness_LS = fitness_LS';
nfes = nfes + popsize_LS;
%%%%%%%%%%%%%

[Sorted_FitVector, Indecis] = sort(fitness_LS);
popLS = popLS(Indecis,:);%sorting the points based on obtaind result
%==========================================================================

%Finding the Best point in the group=======================================
BestPoint = popLS(1, :);
F = Sorted_FitVector(1);%saving the first best fitness
%%%%%%%%%%%%%

run_funcvals = [run_funcvals;fitness];

run_funcvals = [run_funcvals;fitness_LS];



bsf_fit_var = 1e+30;
bsf_index = 0;
bsf_solution = zeros(1, problem_size);

%%%%%%%%%%%%%%%%%%%%%%%% for out
for i = 1 : pop_size
    nfes = nfes + 1;
    
    if fitness(i) < bsf_fit_var
        bsf_fit_var = fitness(i);
        bsf_solution = pop(i, :);
        bsf_index = i;
    end
    
    if nfes > max_nfes; break; end
end
bsf_fit = bsf_fit_var;
bsf_sol = bsf_solution;

%%%%%%%%%%%%%%%%%%%%%%%% for out

memory_sf = 0.5 .* ones(memory_size, 1);
memory_cr = 0.5 .* ones(memory_size, 1);

memory_freq = freq_inti*ones(memory_size, 1);
memory_pos = 1;

archive.NP = round(arc_rate * pop_size); % the maximum size of the archive
archive.pop = zeros(0, problem_size); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions

%% main loop
gg=0;  %%% generation counter used For Sin
igen =1;  %%% generation counter used For LS

flag1 = false;
flag2 = false;

counter = 0;
if isdisplay == 1
    disp(['Generation' , '    Func-count', '   Best-f(x)']);
end
while nfes < max_nfes
    gg=gg+1;
    if isdisplay == 1
        disp(['    ',num2str(gg), '            ',num2str(nfes,'%d'),...
            '        ', num2str(bsf_fit_var,'%04.2f')]);
    end
    pop = popold; % the old population becomes the current population
    [temp_fit, sorted_index] = sort(fitness, 'ascend');
    
    mem_rand_index = ceil(memory_size * rand(pop_size, 1));
    mu_sf = memory_sf(mem_rand_index);
    mu_cr = memory_cr(mem_rand_index);
    mu_freq = memory_freq(mem_rand_index);
    
    %% for generating crossover rate
    cr = normrnd(mu_cr, 0.1);
    term_pos = find(mu_cr == -1);
    cr(term_pos) = 0;
    cr = min(cr, 1);
    cr = max(cr, 0);
    
    %% for generating scaling factor
    sf = mu_sf + 0.1 * tan(pi * (rand(pop_size, 1) - 0.5));
    pos = find(sf <= 0);
    
    while ~ isempty(pos)
        sf(pos) = mu_sf(pos) + 0.1 * tan(pi * (rand(length(pos), 1) - 0.5));
        pos = find(sf <= 0);
    end
    
    
    freq = mu_freq + 0.1 * tan(pi*(rand(pop_size, 1) - 0.5));
    pos_f = find(freq <=0);
    while ~ isempty(pos_f)
        freq(pos_f) = mu_freq(pos_f) + 0.1 * tan(pi * (rand(length(pos_f), 1) - 0.5));
        pos_f = find(freq <= 0);
    end
    
    sf = min(sf, 1);
    freq = min(freq, 1);
    
    if(nfes <= max_nfes/2)
        c=rand;
        if(c<0.5)
            sf = 0.5.*( sin(2.*pi.*freq_inti.*gg+pi) .* ((G_Max-gg)/G_Max) + 1 ) .*...
                ones(pop_size,problem_size);
        else
            sf = 0.5 *( sin(2*pi .* freq(:, ones(1, problem_size)) .* gg) .* ...
                (gg/G_Max) + 1 ) .* ones(pop_size,problem_size);
        end
    end
    
    r0 = [1 : pop_size];
    popAll = [pop; archive.pop];
    fitnessAll = [fitness; archive.funvalues];
    [r1, r2] = gnR1R2(pop_size, size(popAll, 1), r0);
    
    pNP = max(round(p_best_rate * pop_size), 2); %% choose at least two best solutions
    randindex = ceil(rand(1, pop_size) .* pNP); %% select from [1, 2, 3, ..., pNP]
    randindex = max(1, randindex); %% to avoid the problem that rand = 0 and thus ceil(rand) = 0
    pbest = pop(sorted_index(randindex), :); %% randomly choose one of the top 100p% solutions
    
    %% fitness-based F
    id = fitness(r1)-fitnessAll(r2)>max(fitness)-min(fitness);
    [~,id1] = max(fitness);
    tmp = r2;
    tmp(id) = id1;
    sf_r1_r2 = abs(fitness(r1)-fitnessAll(r2))./(max(fitness)-min(fitness(r1)));
    sf_r1_r2 = min(0.999,sf_r1_r2); % value of r1r2
    
    sf_r1_r2_p = sf_r1_r2.^(1/2); % positive and power-1/2
    
    sf_r1_r2_rp = (1-sf_r1_r2).^2; % negative and power-2, better
    
    
    sf_pbest_r3 = abs(fitness(sorted_index(randindex))-fitness(r0))./...
        (max(abs(fitness(sorted_index(randindex))-fitness(r0)))-...
        min(abs(fitness(sorted_index(randindex))-fitness(r0))));
    sf_pbest_r3 = min(0.999,sf_pbest_r3);
    sf_pbest_r3 = max(0.001,sf_pbest_r3); % value of pbest_r3
    
    sf_pbest_r3_p = sf_pbest_r3.^(1/2);  % positive and power-1/2
    
    sf_pbest_r3 = 1 - sf_pbest_r3;
    sf_pbest_r3_rp = (2*sf_pbest_r3-1).^(3);
    sf_pbest_r3_rp = 0.5*sf_pbest_r3_rp + 0.5; % negative and power-3, better
    
    
    sf_r3 = abs(fitness(r0)-min(fitness))./(max(fitness)-min(fitness)); % value of r3
    sf_r3 = 1 - sf_r3;
    sf_r3_rp = (2*sf_r3-1).^(3);
    sf_r3_rp = 0.5*sf_r3_rp + 0.5; % negative and power-3, better
    
    flag_sf = rand(pop_size,1);
    sf3 = (flag_sf-1/3<=0).*sf_r1_r2_rp + (flag_sf-2/3>0).*sf_pbest_r3_rp +...
        ((flag_sf-1/3).*(flag_sf-2/3)<0).*sf_r3_rp;
    %% fitness-based F
    %
    sf = sf(:,1);
    %
    index = sf<sf3;
    
    
    sf(index) = sf(index)*1.1;
    sf = min(sf,1);
    sf = max(sf,0);
    
    index1 = sf>sf3;
    sf(index1) = sf(index1)*0.9;
    sf = min(sf,1);
    sf = max(sf,0);
    
    vi = pop + sf(:, ones(1, problem_size)) .* (pbest - pop + pop(r1, :) - popAll(r2, :));
    vi = boundConstraint(vi, pop, lu);
    
    mask = rand(pop_size, problem_size) > cr(:, ones(1, problem_size)); 
    rows = (1 : pop_size)'; cols = floor(rand(pop_size, 1) * problem_size)+1; 
    jrand = sub2ind([pop_size problem_size], rows, cols); mask(jrand) = false;
    ui = vi; ui(mask) = pop(mask);
    
    children_fitness = zeros(1,pop_size);
    
    if UseParallel == 1
        parfor m = 1:pop_size
            children_fitness(m) =  fx(ui(m,:));
        end
    else
        for m = 1:pop_size
            children_fitness(m) =  fx(ui(m,:));
        end
    end
    
    children_fitness = children_fitness';
    
    
    %%%% To check stagnation
    flag = false;
    bsf_fit_var_old = bsf_fit_var;
    %%%%%%%%%%%%%%%%%%%%%%%% for out
    for i = 1 : pop_size
        nfes = nfes + 1;
        
        if children_fitness(i) < bsf_fit_var
            bsf_fit_var = children_fitness(i);
            bsf_solution = ui(i, :);
            bsf_index = i;
        end
        
        if nfes > max_nfes; break; end
    end
    bsf_fit(gg,1) = bsf_fit_var;
    bsf_sol(gg,1:problem_size) = bsf_solution;
    %%%%%%%%%%%%%%%%%%%%%%%% for out
    
    dif = abs(fitness - children_fitness);
    
    
    %% I == 1: the parent is better; I == 2: the offspring is better
    I = (fitness > children_fitness);
    goodCR = cr(I == 1);
    goodF = sf(I == 1);
    goodFreq = freq(I == 1);
    dif_val = dif(I == 1);
    
    %      isempty(popold(I == 1, :))
    archive = updateArchive(archive, popold(I == 1, :), fitness(I == 1));
    
    [fitness, I] = min([fitness, children_fitness], [], 2);
    
    run_funcvals = [run_funcvals; fitness];
    
    popold = pop;
    popold(I == 2, :) = ui(I == 2, :);
    
    num_success_params = numel(goodCR);
    
    if num_success_params > 0
        sum_dif = sum(dif_val);
        dif_val = dif_val / sum_dif;
        
        %% for updating the memory of scaling factor
        memory_sf(memory_pos) = (dif_val' * (goodF .^ 2)) / (dif_val' * goodF);
        
        %% for updating the memory of crossover rate
        if max(goodCR) == 0 || memory_cr(memory_pos)  == -1
            memory_cr(memory_pos)  = -1;
        else
            memory_cr(memory_pos) = (dif_val' * (goodCR .^ 2)) / (dif_val' * goodCR);
        end
        
        %% for updating the memory of freq
        if max(goodFreq) == 0 || memory_freq(memory_pos)  == -1
            memory_freq(memory_pos)  = -1;
        else
            memory_freq(memory_pos) = (dif_val' * (goodFreq .^ 2)) / (dif_val' * goodFreq);
        end
        
        memory_pos = memory_pos + 1;
        if memory_pos > memory_size;  memory_pos = 1; end
    end
    
    %% for resizing the population size
    % It is supposed the solution is nearly globelly optimal now,
    % This part is to decrese the population size to reudce computational cost.
    if (nfes >= (max_nfes/2)) && (isresize == 1)
        counter = counter + 1;
        plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
        
        if pop_size > plan_pop_size
            reduction_ind_num = pop_size - plan_pop_size;
            if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
            
            if counter == 1
                count = 0;
                stop = false;
                niches_no = 0;
                %% Change here Niche-based reduction
                %% Firstly exculde the best niche (Half of the individuals)
                %% Step 1: sort according fitness to pick up the best individual
                [valBest indBest] = sort(fitness, 'ascend');  %%descend
                best_ind = indBest(1);
                best_mem = pop(best_ind,:);
                
                %% Step 2: find E-distance between best_mem and all others
                %% To Choose neighbourhood region to the best individual
                Dis = pdist2(pop,best_mem,'euclidean'); % euclidean distance
                
                %% Sort and chooose smallest distance to have higher diversity
                [Dis_ordered idx_ordered] = sort(Dis, 'ascend');
                best_niche_size = round(pop_size/2);
                
                %% Select the memebrs of the best niche
                best_niche = pop(idx_ordered(1:best_niche_size), :); %%% including best also so start from 1
                best_niche_idx = idx_ordered(1:best_niche_size);
                best_niche_fit = fitness(idx_ordered(1:best_niche_size));
                
                %% Delete them temporaily
                pop(idx_ordered(1:best_niche_size), :) = [];
                popold(idx_ordered(1:best_niche_size), :) = [];
                fitness(idx_ordered(1:best_niche_size)) = [];
                
                niche_size = 20;
                
                %% Define temp arrays to store the remainning individuals
                pop_temp = [];
                popold_temp = [];
                fitness_temp = [];
                
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');  %%descend
                    best_ind = indBest(1);
                    best_mem = pop(best_ind,:);
                    %                   fitnesshere = fitness(best_ind)
                    
                    Dis = pdist2(pop,best_mem,'euclidean'); % euclidean distance
                    %% Sort and chooose smallest distance to have higher diversity
                    [Dis_ordered idx_ordered] = sort(Dis, 'ascend');
                    
                    [NP_curr, D] = size(pop);
                    if(NP_curr < niche_size)
                        niche_size = NP_curr;
                    end
                    
                    niche = pop(idx_ordered(1:niche_size), :); %%% including best also so start from 1
                    niche_fitness = fitness(idx_ordered(1:niche_size), :);
                    niche_idx = idx_ordered(1:niche_size);
                    
                    niches_no = niches_no + 1;
                    
                    niche_idx_keep = [];
                    niche_idx_delete = [];
                    
                    %% Now remove half of them exculding best
                    %% The best way to remove is to mark those individuals by -1
                    for t = 2 : ((niche_size/2)+1)  %% To not include best and remove it, start from 2 not 1!
                        %                       fitness(idx_ordered(t)) = -1;
                        niche_idx_delete = [niche_idx_delete; t];
                        
                        count = count + 1;
                        if count == reduction_ind_num
                            stop = true;
                            break;
                        end
                    end
                    
                    %% Then those indecies should be removed from pop temporaily
                    niche(niche_idx_delete, :) = [];
                    niche_fitness(niche_idx_delete, :) = [];
                    
                    %% Keep them to not get lost
                    pop_temp = [pop_temp; niche];
                    popold_temp = [popold_temp; niche];
                    fitness_temp = [fitness_temp; niche_fitness];
                    
                    %% And then remove from pop to start with new niche
                    pop(niche_idx_delete, :) = [];
                    popold(niche_idx_delete, :) = [];
                    fitness(niche_idx_delete, :) = [];
                    
                    if stop == true
                        break;
                    end
                end
                
                %% Check again as one element is -1!!
                pos_fit = find(fitness == -1);
                while ~ isempty(pos_fit)
                    fitness(pos_fit) = [];
                    popold(pos_fit,:) = [];
                    pop(pos_fit,:) = [];
                    pos_fit = find(fitness == -1);
                end
                
                pop = [pop; best_niche];
                popold = [popold; best_niche];
                fitness = [fitness; best_niche_fit];
                
                pop_size = pop_size - reduction_ind_num;
            else
                
                pop_size = pop_size - reduction_ind_num;
                SEL = round(ps*pop_size);
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind,:) = [];
                    pop(worst_ind,:) = [];
                    fitness(worst_ind,:) = [];
                end
            end
            archive.NP = round(arc_rate * pop_size);
            
            if size(archive.pop, 1) > archive.NP
                rndpos = randperm(size(archive.pop, 1));
                rndpos = rndpos(1 : archive.NP);
                archive.pop = archive.pop(rndpos, :);
                archive.funvalues = archive.funvalues(rndpos, :);
            end
        end
        %           pop_size = pop_size
    end
    
    %%%%%%%%%%%%%%% Call LS based on Gaussian works when NP is less than 20 for the first time  %%%%%
    if pop_size <= 20
        counter = counter + 1;
    end
    
    if counter == 1
        flag_LS = true;
    else
        flag_LS = false;
    end
    
    if flag_LS == true
        r_index = randi([1 pop_size],1,popsize_LS);
        %%% Pick 10 random individuals from L-SHADE pop
        for gen_LS = 0 : GenMaxSelected
            New_Point = [];%creating new point
            FitVector = [];%creating vector of fitness functions
            
            for i = 1 : popsize_LS
                [NP, fit] = LS_Process(fx,popLS(i,:),S,igen,BestPoint,UseParallel);
                New_Point = [New_Point;NP];
                FitVector = [FitVector,fit];
            end
            
            %%%%
            fittemp = FitVector;
            for i = 1 : popsize_LS
                %%% Update those 10 random individuals from pop L-SHADE
                if FitVector(i) < fitness(r_index(i))
                    fitness (r_index(i)) = FitVector(i);
                    pop(r_index(i),:) = New_Point(i,:);
                    
                else
                    fittemp(i) =  fitness (r_index(i));
                end
                
                %%%% Update best individual L-SHADE
                if FitVector(i) < bsf_fit_var
                    bsf_fit_var = FitVector(i);
                    bsf_solution = New_Point(i,:);
                end
                
                nfes = nfes + 1;
                if nfes > max_nfes; break; end
            end
            
            bsf_fit(gg,1) = bsf_fit_var;
            bsf_sol(gg,1:problem_size) = bsf_solution;
            
            %%%%%% To recored those changes
            fittemp = fittemp';
            run_funcvals = [run_funcvals; fittemp];
            
            %%%%%%%%%%%%%%
            
            [SortedFit,SortedIndex] = sort(FitVector);
            New_Point = New_Point(SortedIndex,:);
            BestPoint = New_Point(1,:);%first point is the best
            BestFitness = SortedFit(1,1);
            popLS = New_Point;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%%%%nfes

end

function archive = updateArchive(archive, pop, funvalue)
% Update the archive with input solutions
%   Step 1: Add new solution to the archive
%   Step 2: Remove duplicate elements
%   Step 3: If necessary, randomly remove some solutions to maintain the archive size
%
% Version: 1.1   Date: 2008/04/02
% Written by Jingqiao Zhang (jingqiao@gmail.com)

if archive.NP == 0, return; end

if size(pop, 1) ~= size(funvalue,1), error('check it'); end

% Method 2: Remove duplicate elements
popAll = [archive.pop; pop ];
funvalues = [archive.funvalues; funvalue ];
[dummy IX]= unique(popAll, 'rows');
if length(IX) < size(popAll, 1) % There exist some duplicate solutions
  popAll = popAll(IX, :);
  funvalues = funvalues(IX, :);
end

if size(popAll, 1) <= archive.NP   % add all new individuals
  archive.pop = popAll;
  archive.funvalues = funvalues;
else                % randomly remove some solutions
  rndpos = randperm(size(popAll, 1)); % equivelent to "randperm";
  rndpos = rndpos(1 : archive.NP);
  
  archive.pop = popAll  (rndpos, :);
  archive.funvalues = funvalues(rndpos, :);
end

end

%This function is used as a local search, and creates some
%new points based on Gaussian Walks.

%**************************************************************************
%The input function is:                                                   %
%Point: the input point which is going to be diffused                     %
%S: structure of problem information                                      %
%g: generation number                                                     %
%BestPoint: the best point in group                                       %                               %
%==========================================================================
%The output function is:                                                  %
%createPoint: the new points created by Diffusion process                 %
%fitness: the value of fitness function                                   %
%**************************************************************************

function [createPoint, fitness] = LS_Process(fx,Point,S,g,BestPoint,UseParallel)
   

    GeneratePoint = normrnd(BestPoint, (log(g)/g)*(abs((Point - BestPoint))), [1 size(Point,2)]) + ...
        (randn*BestPoint - randn*Point);
    
    %check bounds of generated point
    GeneratePoint = Bound_Checking(GeneratePoint,S.Lband,S.Uband);
    
%     size(GeneratePoint)  
%     for i=1:size(Point,2)
%         if GeneratePoint(1,i) > S.Uband
%             fprintf('violate upper');
%         end
%         if GeneratePoint(1,i) < S.Lband
%              fprintf('violate lower');
%         end  
%     end
    
%     fitness = feval(fhd,GeneratePoint',S.FuncNo);
    fitness = zeros(1,length(GeneratePoint(:,1)));
    
    if UseParallel == 1
        parfor m = 1:length(GeneratePoint(:,1))
            fitness(m) =  fx(GeneratePoint(m,:));
        end
    else
        for m = 1:length(GeneratePoint(:,1))
            fitness(m) =  fx(GeneratePoint(m,:));
        end
    end

    createPoint = GeneratePoint;
    %======================================================================
end

function G_Max = G_Max_calc(problem_size,max_nfes)

lb = -100*ones(1,problem_size); % 下界，行向量
ub = 100*ones(1,problem_size); % 上界，行向量
lu = [lb; ub];
G_Max0 = max_nfes; % 根据目标函数计算次数，随时调整
%% 问题参数

%% DE本身参数
freq_inti = 0.5;
% max_nfes = 30000;%5000 * problem_size; % 目标函数计算次数
pb = 0.4;
ps = .5;
S.Ndim = problem_size;
S.Lband = lb;
S.Uband = ub;
GenMaxSelected = 250; %%% For local search
%% DE本身参数


%% parameter settings for L-SHADE
p_best_rate = 0.11;    %0.11
arc_rate = 1.4;
memory_size = 5;
pop_size = 18 * problem_size;   %18*D
SEL = round(ps*pop_size);

max_pop_size = pop_size;
min_pop_size = 4.0;

run_funcvals = [];

nfes = 0;
%% Initialize the main population
popold = repmat(lu(1, :), pop_size, 1) + rand(pop_size, problem_size) .*...
    (repmat(lu(2, :) - lu(1, :), pop_size, 1));
pop = popold; % the old population becomes the current population

fitness = zeros(1,pop_size);
for m = 1:pop_size
    
    fitness(m) = 1;
    
end
fitness = fitness';

%%% Initialize LS population
flag_LS = false;
counter = 0;
popsize_LS = 10;

%%% Initialize LS population to re-start them
popLS = repmat(lu(1, :), popsize_LS, 1) + rand(popsize_LS, problem_size) .* ...
    (repmat(lu(2, :) - lu(1, :), popsize_LS, 1));
fitness_LS = zeros(1,popsize_LS);
for m = 1:popsize_LS
    
    fitness_LS(m) = 1;
    
end

fitness_LS = fitness_LS';
nfes = nfes + popsize_LS;
%%%%%%%%%%%%%

[Sorted_FitVector, Indecis] = sort(fitness_LS);
popLS = popLS(Indecis,:);%sorting the points based on obtaind result
%==========================================================================

%Finding the Best point in the group=======================================
BestPoint = popLS(1, :);
F = Sorted_FitVector(1);%saving the first best fitness
%%%%%%%%%%%%%

run_funcvals = [run_funcvals;fitness];

run_funcvals = [run_funcvals;fitness_LS];



bsf_fit_var = 1e+30;
bsf_index = 0;
bsf_solution = zeros(1, problem_size);

%%%%%%%%%%%%%%%%%%%%%%%% for out
for i = 1 : pop_size
    nfes = nfes + 1;
    
    if fitness(i) < bsf_fit_var
        bsf_fit_var = fitness(i);
        bsf_solution = pop(i, :);
        bsf_index = i;
    end
    
    if nfes > max_nfes; break; end
end
%%%%%%%%%%%%%%%%%%%%%%%% for out

memory_sf = 0.5 .* ones(memory_size, 1);
memory_cr = 0.5 .* ones(memory_size, 1);

memory_freq = freq_inti*ones(memory_size, 1);
memory_pos = 1;

archive.NP = round(arc_rate * pop_size); % the maximum size of the archive
archive.pop = zeros(0, problem_size); % the solutions stored in te archive
archive.funvalues = zeros(0, 1); % the function value of the archived solutions

%% main loop
gg=0;  %%% generation counter used For Sin
igen =1;  %%% generation counter used For LS

flag1 = false;
flag2 = false;

counter = 0;
while nfes < max_nfes
    gg=gg+1;
%     [gg, bsf_fit_var]
    pop = popold; % the old population becomes the current population
    ui = popold; 
    
    children_fitness = zeros(1,pop_size);
    for m = 1:pop_size
        
        children_fitness(m) =  1;
        
    end
    children_fitness = children_fitness';
    
    
    %%%% To check stagnation
    flag = false;
    bsf_fit_var_old = bsf_fit_var;
    %%%%%%%%%%%%%%%%%%%%%%%% for out
    for i = 1 : pop_size
        nfes = nfes + 1;
        
        if children_fitness(i) < bsf_fit_var
            bsf_fit_var = children_fitness(i);
            bsf_solution = ui(i, :);
            bsf_index = i;
        end
        
        if nfes > max_nfes; break; end
    end
    %%%%%%%%%%%%%%%%%%%%%%%% for out
    
    %% for resizing the population size
    if (nfes >= (max_nfes/2))
        counter = counter + 1;
        plan_pop_size = round((((min_pop_size - max_pop_size) / max_nfes) * nfes) + max_pop_size);
        
        if pop_size > plan_pop_size
            reduction_ind_num = pop_size - plan_pop_size;
            if pop_size - reduction_ind_num <  min_pop_size; reduction_ind_num = pop_size - min_pop_size;end
            
            if counter == 1
                count = 0;
                stop = false;
                niches_no = 0;
                %% Change here Niche-based reduction
                %% Firstly exculde the best niche (Half of the individuals)
                %% Step 1: sort according fitness to pick up the best individual
                [valBest indBest] = sort(fitness, 'ascend');  %%descend
                best_ind = indBest(1);
                best_mem = pop(best_ind,:);
                
                %% Step 2: find E-distance between best_mem and all others
                %% To Choose neighbourhood region to the best individual
                Dis = pdist2(pop,best_mem,'euclidean'); % euclidean distance
                
                %% Sort and chooose smallest distance to have higher diversity
                [Dis_ordered idx_ordered] = sort(Dis, 'ascend');
                best_niche_size = round(pop_size/2);
                
                %% Select the memebrs of the best niche
                best_niche = pop(idx_ordered(1:best_niche_size), :); %%% including best also so start from 1
                best_niche_idx = idx_ordered(1:best_niche_size);
                best_niche_fit = fitness(idx_ordered(1:best_niche_size));
                
                %% Delete them temporaily
                pop(idx_ordered(1:best_niche_size), :) = [];
                popold(idx_ordered(1:best_niche_size), :) = [];
                fitness(idx_ordered(1:best_niche_size)) = [];
                
                niche_size = 20;
                
                %% Define temp arrays to store the remainning individuals
                pop_temp = [];
                popold_temp = [];
                fitness_temp = [];
                
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');  %%descend
                    best_ind = indBest(1);
                    best_mem = pop(best_ind,:);
                    %                   fitnesshere = fitness(best_ind)
                    
                    Dis = pdist2(pop,best_mem,'euclidean'); % euclidean distance
                    %% Sort and chooose smallest distance to have higher diversity
                    [Dis_ordered idx_ordered] = sort(Dis, 'ascend');
                    
                    [NP_curr, D] = size(pop);
                    if(NP_curr < niche_size)
                        niche_size = NP_curr;
                    end
                    
                    niche = pop(idx_ordered(1:niche_size), :); %%% including best also so start from 1
                    niche_fitness = fitness(idx_ordered(1:niche_size), :);
                    niche_idx = idx_ordered(1:niche_size);
                    
                    niches_no = niches_no + 1;
                    
                    niche_idx_keep = [];
                    niche_idx_delete = [];
                    
                    %% Now remove half of them exculding best
                    %% The best way to remove is to mark those individuals by -1
                    for t = 2 : ((niche_size/2)+1)  %% To not include best and remove it, start from 2 not 1!
                        %                       fitness(idx_ordered(t)) = -1;
                        niche_idx_delete = [niche_idx_delete; t];
                        
                        count = count + 1;
                        if count == reduction_ind_num
                            stop = true;
                            break;
                        end
                    end
                    
                    %% Then those indecies should be removed from pop temporaily
                    niche(niche_idx_delete, :) = [];
                    niche_fitness(niche_idx_delete, :) = [];
                    
                    %% Keep them to not get lost
                    pop_temp = [pop_temp; niche];
                    popold_temp = [popold_temp; niche];
                    fitness_temp = [fitness_temp; niche_fitness];
                    
                    %% And then remove from pop to start with new niche
                    pop(niche_idx_delete, :) = [];
                    popold(niche_idx_delete, :) = [];
                    fitness(niche_idx_delete, :) = [];
                    
                    if stop == true
                        break;
                    end
                end
                
                %% Check again as one element is -1!!
                pos_fit = find(fitness == -1);
                while ~ isempty(pos_fit)
                    fitness(pos_fit) = [];
                    popold(pos_fit,:) = [];
                    pop(pos_fit,:) = [];
                    pos_fit = find(fitness == -1);
                end
                
                pop = [pop; best_niche];
                popold = [popold; best_niche];
                fitness = [fitness; best_niche_fit];
                
                pop_size = pop_size - reduction_ind_num;
            else
                
                pop_size = pop_size - reduction_ind_num;
                SEL = round(ps*pop_size);
                for r = 1 : reduction_ind_num
                    [valBest indBest] = sort(fitness, 'ascend');
                    worst_ind = indBest(end);
                    popold(worst_ind,:) = [];
                    pop(worst_ind,:) = [];
                    fitness(worst_ind,:) = [];
                end
            end
            archive.NP = round(arc_rate * pop_size);
            
            if size(archive.pop, 1) > archive.NP
                rndpos = randperm(size(archive.pop, 1));
                rndpos = rndpos(1 : archive.NP);
                archive.pop = archive.pop(rndpos, :);
                archive.funvalues = archive.funvalues(rndpos, :);
            end
        end
        %           pop_size = pop_size
    end
    
    %%%%%%%%%%%%%%% Call LS based on Gaussian works when NP is less than 20 for the first time  %%%%%
    if pop_size <= 20
        counter = counter + 1;
    end
    
    if counter == 1
        flag_LS = true;
    else
        flag_LS = false;
    end
    
    if flag_LS == true
        r_index = randi([1 pop_size],1,popsize_LS);
        %%% Pick 10 random individuals from L-SHADE pop
        for gen_LS = 0 : GenMaxSelected
            New_Point = [];%creating new point
            FitVector = [];%creating vector of fitness functions
            
            for i = 1 : popsize_LS
                
                NP = pop(1,:);
                fit = 1;
                New_Point = [New_Point;NP];
                FitVector = [FitVector,fit];
            end
            
            %%%%
            fittemp = FitVector;
            for i = 1 : popsize_LS
                %%% Update those 10 random individuals from pop L-SHADE
                if FitVector(i) < fitness(r_index(i))
                    fitness (r_index(i)) = FitVector(i);
                    pop(r_index(i),:) = New_Point(i,:);
                    
                else
                    fittemp(i) =  fitness (r_index(i));
                end
                
                %%%% Update best individual L-SHADE
                if FitVector(i) < bsf_fit_var
                    bsf_fit_var = FitVector(i);
                    bsf_solution = New_Point(i,:);
                end
                
                nfes = nfes + 1;
                if nfes > max_nfes; break; end
            end
            
            %%%%%% To recored those changes
            fittemp = fittemp';
            run_funcvals = [run_funcvals; fittemp];
            
            %%%%%%%%%%%%%%
            
            [SortedFit,SortedIndex] = sort(FitVector);
            New_Point = New_Point(SortedIndex,:);
            BestPoint = New_Point(1,:);%first point is the best
            BestFitness = SortedFit(1,1);
            popLS = New_Point;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end %%%%%%%%nfes
G_Max = gg;
end

function [r1, r2] = gnR1R2(NP1, NP2, r0)

% gnA1A2 generate two column vectors r1 and r2 of size NP1 & NP2, respectively
%    r1's elements are choosen from {1, 2, ..., NP1} & r1(i) ~= r0(i)
%    r2's elements are choosen from {1, 2, ..., NP2} & r2(i) ~= r1(i) & r2(i) ~= r0(i)
%
% Call:
%    [r1 r2 ...] = gnA1A2(NP1)   % r0 is set to be (1:NP1)'
%    [r1 r2 ...] = gnA1A2(NP1, r0) % r0 should be of length NP1
%
% Version: 2.1  Date: 2008/07/01
% Written by Jingqiao Zhang (jingqiao@gmail.com)

NP0 = length(r0);

r1 = floor(rand(1, NP0) * NP1) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = (r1 == r0);
    if sum(pos) == 0
        break;
    else % regenerate r1 if it is equal to r0
        r1(pos) = floor(rand(1, sum(pos)) * NP1) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r1 in 1000 iterations');
    end
end

r2 = floor(rand(1, NP0) * NP2) + 1;
%for i = 1 : inf
for i = 1 : 99999999
    pos = ((r2 == r1) | (r2 == r0));
    if sum(pos)==0
        break;
    else % regenerate r2 if it is equal to r0 or r1
        r2(pos) = floor(rand(1, sum(pos)) * NP2) + 1;
    end
    if i > 1000, % this has never happened so far
        error('Can not genrate r2 in 1000 iterations');
    end
end
end

%This function is used for L-SHADE bound checking 
function vi = boundConstraint (vi, pop, lu)

% if the boundary constraint is violated, set the value to be the middle
% of the previous value and the bound
%

[NP, D] = size(pop);  % the population size and the problem's dimension

%% check the lower bound
xl = repmat(lu(1, :), NP, 1);

pos = vi < xl;
vi(pos) = (pop(pos) + xl(pos)) / 2;
% vi(pos) = xl(pos);
%% check the upper bound
xu = repmat(lu(2, :), NP, 1);
pos = vi > xu;
vi(pos) = (pop(pos) + xu(pos)) / 2;
% vi(pos) = xu(pos);
end

%This function is used for LS bound checking 
function p = Bound_Checking(p,lowB,upB)
    for i = 1 : size(p,1)
        upper = double(gt(p(i,:),upB));
        lower = double(lt(p(i,:),lowB));
        up = find(upper == 1);
        lo = find(lower == 1);
        if (size(up,2)+ size(lo,2) > 0 )
            for j = 1 : size(up,2)
%                 fprintf('here');
                p(i, up(j)) = (upB(up(j)) - lowB(up(j)))*rand()...
                    + lowB(up(j));
            end
            for j = 1 : size(lo,2)
                p(i, lo(j)) = (upB(lo(j)) - lowB(lo(j)))*rand()...
                    + lowB(lo(j));
            end
        end
    end
end