%% Talha Akhlaq
% ECE300 Communication Theory
% Problem Set 2: Conditional Distributions & Decision Theory
% Question 4
clc; clear; close all;

RED = 1; BLUE = 2; % indices 

%% Part (d): 

N = 1e5;
cases = { ...
  struct('b1',2,'r1',8,'b2',2,'r2',8,'name','Case I'), ...
  struct('b1',6,'r1',4,'b2',2,'r2',8,'name','Case II'), ...
  struct('b1',2,'r1',8,'b2',6,'r2',4,'name','Case III'), ...
  struct('b1',6,'r1',4,'b2',6,'r2',4,'name','Case IV')};

fprintf('Question 4: \n\n');
for k = 1:numel(cases)
    C = cases{k};
    T = theory_from_urns(C.r1, C.b1, C.r2, C.b2);

    Perr_ML_th  = rule_error(T.prior, T.L, T.rule_ML);
    Perr_MAP_th = rule_error(T.prior, T.L, T.rule_MAP);

    Perr_ML_emp  = simulate_urns(N, C.r1, C.b1, C.r2, C.b2, T.rule_ML);
    Perr_MAP_emp = simulate_urns(N, C.r1, C.b1, C.r2, C.b2, T.rule_MAP);

    fprintf('%s: r1=%d b1=%d | r2=%d b2=%d\n', C.name, C.r1, C.b1, C.r2, C.b2);
    fprintf('  Priors        P(RedI)=%.4f  P(BlueI)=%.4f\n', T.prior(RED), T.prior(BLUE));
    fprintf('  Likelihoods   P(y|i):  RedII=[%.4f %.4f]  BlueII=[%.4f %.4f]\n', ...
        T.L(RED,RED), T.L(RED,BLUE), T.L(BLUE,RED), T.L(BLUE,BLUE));
    fprintf('  Posteriors    P(i|y):  RedII=[%.4f %.4f]  BlueII=[%.4f %.4f]\n', ...
        T.post(RED,RED), T.post(RED,BLUE), T.post(BLUE,RED), T.post(BLUE,BLUE));
    fprintf('  ML  rule      = [%s %s]\n',  pick(T.rule_ML(RED)),  pick(T.rule_ML(BLUE)));
    fprintf('  MAP rule      = [%s %s]\n', pick(T.rule_MAP(RED)), pick(T.rule_MAP(BLUE)));
    fprintf('  Error theory  ML=%.5f  MAP=%.5f\n',  Perr_ML_th,  Perr_MAP_th);
    fprintf('  Error empirical ML=%.5f  MAP=%.5f  N=%d\n\n', Perr_ML_emp, Perr_MAP_emp, N);
end

% Comparison:
% Results match theory
% Case I: MAP better (priors favor Red)
% Case II: MAP better (priors favor Blue)
% Case III: MAP better (priors favor Red)
% Case IV: MAP better (priors favor Blue)
% MAP = ML only if priors are equal or do not change the decision/argmax
% MAP > ML whenever priors are skewed

%% Parts (a-c): 
function T = theory_from_urns(r1,b1,r2,b2) % Priors, likelihoods, posteriors, ML, MAP
RED = 1; BLUE = 2;
prior = [r1/(r1+b1), b1/(r1+b1)];
den = r2 + b2 + 1;
L = zeros(2,2);                              % rows=y, cols=i
L(RED,RED)   = (r2+1)/den;                   % P(RedII | RedI)
L(BLUE,RED)  =  b2   /den;                   % P(BlueII| RedI)
L(RED,BLUE)  =  r2   /den;                   % P(RedII | BlueI)
L(BLUE,BLUE) = (b2+1)/den;                   % P(BlueII| BlueI)
post = zeros(2,2);                           % Bayes posteriors P(i|y)
for y = [RED BLUE]
    w = prior .* L(y,:);
    post(y,:) = w / sum(w);
end
% Part (a): 1x2 array
[~, rule_ML] = max(L,[],2);                  % argmax_i P(y|i)
W = (ones(2,1)*prior).*L;                    % prior weighting for MAP
[~, rule_MAP] = max(W,[],2);                 % argmax_i pi_i P(y|i)
T = struct('prior',prior,'L',L,'post',post,'rule_ML',rule_ML.', 'rule_MAP',rule_MAP.');
end

function Perr = rule_error(prior,L,rule)     % Theoretical error for deterministic rule
RED = 1; BLUE = 2;
P_y_i = [prior(RED)*L(RED,RED),  prior(BLUE)*L(RED,BLUE); ...
         prior(RED)*L(BLUE,RED), prior(BLUE)*L(BLUE,BLUE)];
acc = P_y_i(RED,  rule(RED)) + P_y_i(BLUE, rule(BLUE));
Perr = 1 - acc;
end

function Perr = simulate_urns(N,r1,b1,r2,b2,rule) % Vectorized randi simulator
RED = 1; BLUE = 2;
u1 = randi(r1+b1, N, 1);                      % Urn I draw
X  = ones(N,1); X(u1>r1) = BLUE;              % 1=RedI, 2=BlueI
den = r2 + b2 + 1;                            % Urn II size after transfer
th  = r2 + (X==RED);                          % red threshold after transfer
u2  = randi(den, N, 1);                       % Urn II draw
Y   = ones(N,1); Y(u2>th) = BLUE;             % 1=RedII, 2=BlueII
hatX = rule(Y).';                             % decision
Perr = mean(hatX ~= X);
end

function s = pick(idx)                        % class name
if idx==1, s='RedI'; else, s='BlueI'; end
end

% Vectorized randi simulates all N trials at once, achieving time
% efficiency (avoids loops)


