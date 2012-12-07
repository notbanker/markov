classdef dsStructUtil
    % Some utilities to translate between datasets and structs and to perform 
    % some operations on structs in the same way they can be done with
    % datasets.
    %
    % NOTE: this DOES index the string fields and is best used when string
    % data is very repetitive or fast lookups on the string data is
    % important.
    
    methods (Static) % Translating between structs and datasets
        
        % Note that if you set ds_ = struct2Dataset(dataset2Struct(ds)) you
        % may have isequal(ds_,ds) = false.  This is because there are
        % NaNs in the dataset.  Try isequalwithequalnans(ds_,ds) instead.
        % The fields may also have been returned in a different order,
        % which will still cause isequalwithequalnans(ds_,ds) to be false.
        
        function st = dataset2Struct(ds)
            % Peel off the variables in a dataset and return them in a
            % struct.  
            %
            % For any string data, create a list of unique string values in the 
            % substruct stringData and return a numerical index in the
            % original field.
            %
            % E.g. if each element of ds.cusip is either '38144LAB6' or
            % '38144LAC4', each element of st.cusip will be either 1 or 2
            % and st.stringData.cusip = {'38144LAB6'; '38144LAC4'}
            dsFields = fieldnames(ds);
            dsFields = dsFields(1:end-1); % the last field is 'Properties'
            st = [];
            
            for f=1:length(dsFields)
                if ischar(ds.(dsFields{f})) || iscell(ds.(dsFields{f}))
                    [st.stringData.(dsFields{f}), ~, st.(dsFields{f})] = unique(ds.(dsFields{f}));                    
                else
                    st.(dsFields{f}) = ds.(dsFields{f});
                end
            end            
        end
        
        function ds = struct2Dataset(st)
            % Convert a struct (possibly with the above stringData substruct) to a dataset.
            if ismember('stringData',fieldnames(st))
                stringData = st.stringData;                
                st = rmfield(st,'stringData');
            else
                stringData = struct;
            end
            
            stFields = fieldnames(st);
            isStringData = ismember(stFields, fieldnames(stringData));
            ds = dataset;
            
            for f=1:length(stFields)
                if isStringData(f)
                    fString = stringData.(stFields{f});
                    ds.(stFields{f}) = fString(st.(stFields{f}));
                else
                    ds.(stFields{f}) = st.(stFields{f});
                end
            end            
        end

        %TODO: do this faster (without using ds as intermediary)
        function ss = struct2SimpleStruct(st)
            warning('struct2SimpleStruct is not optimized for speed.');
            ds = dsStructUtil.struct2Dataset(st);
            ss = dsSimpleStructUtil.dataset2SimpleStruct(ds);
        end

        %TODO: do this faster (without using ds as intermediary)
        function st = simpleStruct2Struct(ss)
            warning('simpleStruct2Struct is not optimized for speed.');
            ds = dsSimpleStructUtil.simpleStruct2Dataset(ss);
            st = dsStructUtil.dataset2Struct(ds);
        end

    end
    
    methods (Static) % Get and set all columns of a dataset for a given set of rows (like datasets can with subsref and subsasgn)
        
        function subSt = getStructSubset(st,rows)   
            % Note: this will return a struct whose stringData field is
            % identical to the original struct, even if fewer strings appear.
            if ismember('stringData',fieldnames(st))
                subSt.stringData = st.stringData;             
                st = rmfield(st,'stringData');
            else
                subSt = [];
            end
            stFields = fieldnames(st);
                
            for f=1:length(stFields)
                subSt.(stFields{f}) = st.(stFields{f})(rows);
            end     
        end
        
        function st = setStructSubset(st,rows,subSt,identicalStringData)
            % If identicaltringData = true, we will assume the stringData fields
            % of st and subSt are identical.  Otherwise, we will check this
            % first and merge them if necessary.
            %
            % Note 1. we assume st and subSt have identical fieldnames.
            % Note 2. identicalStringData = true by default
            if nargin < 4
                identicalStringData = true;
            end
            % Do they have string data?
            if ismember('stringData',fieldnames(st))
                if identicalStringData
                    stringData = st.stringData;
                else    
                    stringFields = fieldnames(st.stringData);
                    for f=1:length(stringFields)
                        % merge the stringData fields and properly handle the
                        % corresponding numeric fields
                        field = st.stringData.(stringFields{f});
                        fieldSub = subSt.stringData.(stringFields{f});
                        mergedField = union(field,fieldSub);  
                        mergedField = mergedField(:); % make sure it's a column 
                        stringData.(stringFields{f}) = mergedField;
                        [~,piFull] = ismember(field, mergedField);  % piFull lists the order that the strings in st's field appear in the merged stringData field
                        [~,piSub] = ismember(fieldSub, mergedField);  
                        
                        % permute the numeric keys as needed
                        nonNaN = ~isnan(st.(stringFields{f}));
                        nonNaNSub = ~isnan(subSt.(stringFields{f}));
                        st.(stringFields{f})(nonNaN) = piFull( st.(stringFields{f})(nonNaN) );
                        subSt.(stringFields{f})(nonNaNSub) = piSub( subSt.(stringFields{f})(nonNaNSub) );
                    end                    
                end
                st.stringData = stringData;                
                stFields = setdiff(fieldnames(st), 'stringData');
            else
                stFields = fieldnames(st);
            end
            
            % insert the appropriate rows
            for f=1:length(stFields)
                st.(stFields{f})(rows) = subSt.(stFields{f});
            end     
        end        
        
        function mergedSt = vertcat(st1,st2,identicalStringData)
            % If identicalStringData = true, we will assume the stringData fields
            % of st1 and st2 are identical.  Otherwise, we will check this
            % first and merge them if necessary.
            %
            % Note 1. we assume st1 and st2 have identical fieldnames.
            % Note 2. identicalStringData = false by default
            if nargin < 3
                identicalStringData = false;
            end
            % Handle some trivial cases first
            if isempty(st1)
                mergedSt = st2;
            elseif isempty(st2)
                mergedSt = st1;
            else
                st1Length = dsStructUtil.length(st1);
                st2Length = dsStructUtil.length(st2);
                mergedSt = dsStructUtil.addNaNs(st1, st2Length);
                mergedSt = dsStructUtil.setStructSubset(mergedSt, st1Length+1:st1Length+st2Length, st2, identicalStringData);
            end
        end            
            
        
    end
    
    methods (Static) % simple utilities for dataset-like structs
        
        function L = length(st)
            % Returns the length of the dataset-like struct
            stFields = setdiff(fieldnames(st),'stringData');
            L = length(st.(stFields{1}));
        end
        
        function st = addNaNs(st,numRows)
            % Adds numRows many NaNs to each field of st
            stFields = setdiff(fieldnames(st),'stringData');
            for f=1:length(stFields)
                st.(stFields{f}) = [st.(stFields{f}); nan(numRows,1)];
            end
        end
        
        function st = cleanNaNs(st)
            % Removes rows in which every entry is NaN (if only some of
            % them are NaN, it leaves the row intact)            
            stFields = setdiff(fieldnames(st),'stringData');
            L = dsStructUtil.length(st);
            F = length(stFields);
            nanFlag = false(L,F);
            
            % First loop through each field and flag the NaNs
            for f=1:F
                nanFlag(:,f) = isnan(st.(stFields{f}));
            end
            rowsToClean = all(nanFlag,2);
            
            % Now loop through each field and dlete the specified rows
            for f=1:F
                st.(stFields{f}) = st.(stFields{f})(~rowsToClean);
            end            
        end
        
        function st = sanitizeStruct(st)
            % Replace any NaNs in fields that came from strings with ones,
            % so struct2Dataset doesn't die a horrible death.
            if isfield(st,'stringData')
                stringFields = fieldnames(st.stringData);
                for f=1:length(stringFields)
                    badEntries = isnan(st.(stringFields{f}));
                    st.(stringFields{f})(badEntries) = ones(nnz(badEntries),1);
                end
            end       
        end
        
        function st = sortrows(st,fieldToSortBy,howToSort)
            % Sort the struct by the values in fieldToSortBy.
            % howToSort = 'ascend' (default) or 'descend'
            if nargin < 3
                howToSort = 'ascend';
            end
            [~,pi_] = sort(st.(fieldToSortBy), howToSort);
            
            stFields = setdiff(fieldnames(st),'stringData');
            for f=1:length(stFields)
                st.(stFields{f}) = st.(stFields{f})(pi_);
            end
        end
        
    end        
    
end

