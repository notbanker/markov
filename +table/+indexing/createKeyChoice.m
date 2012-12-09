function [tildeK,tildeU,tildeR,phiF] = createKeyChoice(f,varargin)
% function [tildeK,tildeU,tildeR,phiF] = createKeyChoice(f,varargin)
%
% Create unique keys into a collection of equal length
% variables (treated as table columns), and lift the user
% supplied mapping f to the level of these keys.
%
% For a more detailed explanation see the docs
%
% tildeK      - Unique key for rows of X
% tildeU(i,:) - The i'th unique row-combinations of column primary keys
% tildeR{i}   - Mapping from columwise primary keys into unique representatives of X(:,i) and f(X)(:,i)
% phiF(k)     - The key of f(X(i,:)) when X(i,:) has key k

% Step 1: Create unique key (i.e. lookup) into rows of X
[K,U,R,O] = table.indexing.createKey(varargin{:});
J = size(U,1);

% Step 2: Apply user choice function f and re-index
n = length(varargin);
allVars = cell(n,1);
[allVars{:}] = applyAndAppend(f,U,R,O,varargin{:});
[hatK,tildeU,tildeR,hatO] = table.indexing.createKey(allVars{:});
tildeK = hatK(K);

debuggingOn=false;
if debuggingOn,
    % Verify that old and new indexes agree
    m = length(varargin{1});
    for m_=1:m,
        for n_=1:n,
            if isnumeric(tildeR{n_}) || islogical(tildeR{n_}),
                vTilde =  tildeR{n_}(tildeU(tildeK(m_),n_));
                v =  R{n_}(U(K(m_),n_));
                assert(isequal(v,vTilde));
            else
                vTilde =  tildeR{n_}{tildeU(tildeK(m_),n_)};
                v =  R{n_}{U(K(m_),n_)};
                assert(isequal(v,vTilde));
            end
        end
    end
    % Verify definitions
    table.assertRightInverse(K,O,(1:length(K)));
    table.assertRightInverse(hatK,hatO,(1:length(hatK)));
end

% Step 3: Lift the user choice function to the new keys
F = [((J+1):2*J)';(1:J)'];
invHatK = table.indexing.rightInverse(hatK,1:J);
phiF = table.indexing.compose(hatK,F,invHatK);

debuggingOn=false;
if debuggingOn,
    % This verifies the commutative diagram (see notes)
    assert(isequalwithequalnans(hatK(F(1:J)),phiF(hatK(1:J))));
    invHatKCheck = table.rightInverseSimple(hatK(1:J));
    assert(isequalwithequalnans(invHatKCheck,invHatK));
    table.assertRightInverse(hatK(1:J),invHatK,(1:J));
end


    function varargout = applyAndAppend(f,U,R,O,varargin)
        
        % Preallocate varargout
        n = length(varargin);
        nJ = length(O);
        nUniqueKeys = length(U);
        for k=1:n,
            if isnumeric(varargin{k}) || islogical(varargin{k}),
                varargout{k} = [varargin{k}(O);zeros(nJ,1)];
            else
                varargout{k} = [varargin{k}(O);cell(nJ,1)];
            end
        end
        
        % In the second half, overwrite the value for the first variable
        for j=1:nUniqueKeys,
            vars = cell(n,1);
            for vNo=1:n,
                if isnumeric(R{vNo}) || islogical(R{vNo}),
                    vars{vNo} = R{vNo}(U(j,vNo));
                else
                    vars{vNo} = R{vNo}{U(j,vNo)};
                end
            end
            varsout_ = cell(n,1);
            [varsout_{:}] = feval(f,vars{:});
            for k=1:n,
                if isnumeric(varargout{k}) || islogical(varargout{k}),
                    varargout{k}(nJ+j) = varsout_{k};
                else
                    varargout{k}{nJ+j} = varsout_{k};
                end
            end
        end
    end
end
