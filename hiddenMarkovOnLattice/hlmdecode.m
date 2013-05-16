function [x,p,lattice] = hlmdecode(y,w,v,varargin) 
% HLMDECODE - Infer posterior state estimate from observations, samples
% of state changes and samples of measurement errors
% 
%   X = HLMDECODE(Y,W,V) uses samples W of changes in a hidden state x
%   and samples of a measurement error V to estimate the mean of the Bayesian posterior
%   hidden state for a very simple discrete time model:
%
%   x(k+1) = x(k) + W
%   y(k)   = x(k) + V
%
%   [X,P,lattice] = HLMDECODE(Y,W,V) also returns the Bayesian posterior
%   probabilities for x taking each of the values in LATTICE. The lattice
%   comprises 2m+1 points and if X is length n P is n x 2m+1 (?)
%
%   X = HLMDECODE(Y,W,V,LATTICE) specifies the lattice points to use. The
%   vector LATTICE should have ascending, evenly spaced entries. If lattice
%   is an integer m then a lattice of size 2m+1 will be created


%% Create HMM model
[lattice,m] = hlmlattice(y,varargin{:});
centeredLattice = lattice-lattice(m+1);
w = w-mean(w);
v = v-mean(v);
W = discretizer.atomic(centeredLattice,w); % weights on lattice representing changes in state
V = discretizer.atomic(centeredLattice,v); % weights on lattice representing measurement errors
tr = translate(W,m+1);
em = translate(V,m+1);
[yNearest,seq] = discretizer.nearest(y,lattice);

%% Decode it
p = hmmdecode(seq',tr,em);
x = lattice*p;


    function em = translate(em0,s0)
        minimumEmitProbability = 1e-8;
        
        % Generate translation invariant hidden Markov model emissions (or transitions)
        n = length(em0);
        m = (n-1)/2;
        em = zeros(n,n);
        for s1=1:n,
            em(s1,:) = shift(em0,s0,s1)';
        end
        
        em = em+minimumEmitProbability;
        em = discretizer.normalizeRows(em);
        
        function em1 = shift(em0,s0,s1)
            % em0(j) is the emission probability to signal j given that we are in state s0
            % em1(j) is the emission probability to signal j given that we are in state s1
            if s1==s0,
                em1 = em0;
            elseif s0>s1,
                n = length(em0);
                em1 = zeros(n,1);
                em1(1:n-(s0-s1)) = em0(1+s0-s1:n);
            elseif s1>s0,
                n = length(em0);
                em1 = zeros(n,1);
                em1(s1-s0+1:end) = em0(1:n-(s1-s0));
            end
        end
    end

end