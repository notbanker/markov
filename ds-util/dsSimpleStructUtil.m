classdef dsSimpleStructUtil
    % Some utilities to translate between datasets and structs and to perform 
    % some operations on structs in the same way they can be done with
    % datasets.
    %
    % NOTE: this DOES NOT index the string fields and is best used for
    % construction, in cases where the string data is not repetitive or
    % when few lookups are done on the string data.
    
    methods (Static) % Translating between structs and datasets
        
        % Note that if you set ds_ = simpleStruct2Dataset(dataset2SimpleStruct(ds)) you
        % may have isequal(ds_,ds) = false.  This is because there are
        % NaNs in the dataset.  Try isequalwithequalnans(ds_,ds) instead.
        % The fields may also have been returned in a different order,
        % which will still cause isequalwithequalnans(ds_,ds) to be false.
        
        function ss = dataset2SimpleStruct(ds)
            % Peel off the variables in a dataset and return them in a
            % struct.  
            dsFields = fieldnames(ds);
            dsFields = dsFields(1:end-1); % the last field is 'Properties'
            ss = [];
            
            for f=1:length(dsFields)
                ss.(dsFields{f}) = ds.(dsFields{f});
            end            
        end
        
        function ds = simpleStruct2Dataset(ss)
            % Convert a struct to a dataset.
            ssFields = fieldnames(ss);
            ds = dataset;
            
            for f=1:length(ssFields)
                ds.(ssFields{f}) = ss.(ssFields{f});
            end            
        end

        function ss = struct2simpleStruct(st)
            ss = dsStructUtil.struct2simpleStruct(st);
        end

        function st = simpleStruct2struct(ss)
            st = dsStructUtil.simpleStruct2struct(ss);
        end

    end

    methods (Static) % Get and set all columns of a dataset for a given set of rows (like datasets can with subsref and subsasgn)
        
        function subSs = getSimpleStructSubset(ss,rows)
            ssFields = fieldnames(ss);
            subSs = [];
                
            for f=1:length(ssFields)
                subSs.(ssFields{f}) = ss.(ssFields{f})(rows);
            end     
        end
        
        function ss = setSimpleStructSubset(ss,rows,subSs)
            % Note: we assume ss and subSs have identical fieldnames.
            ssFields = fieldnames(ss);

            % insert the appropriate rows
            for f=1:length(ssFields)
                ss.(ssFields{f})(rows) = subSs.(ssFields{f});
            end
        end

        function mergedSs = vertcat(ss1,ss2)
            % Note: we assume ss1 and ss2 have identical fieldnames.

            % Handle some trivial cases first
            if isempty(ss1)
                mergedSs = ss2;
            elseif isempty(ss2)
                mergedSs = ss1;
            else
                ss1Length = dsSimpleStructUtil.length(ss1);
                ss2Length = dsSimpleStructUtil.length(ss2);
                mergedSs = dsSimpleStructUtil.addNaNs(ss1, ss2Length);
                mergedSs = dsSimpleStructUtil.setSimpleStructSubset(mergedSs, ss1Length+1:ss1Length+ss2Length, ss2);
            end
        end            
            
        
    end
    
    methods (Static) % simple utilities for dataset-like structs
        
        function L = length(ss)
            % Returns the length of the dataset-like struct
            ssFields = fieldnames(ss);
            L = length(ss.(ssFields{1}));
        end
        
        function ss = addNaNs(ss,numRows)
            % Adds numRows many NaNs to each field of ss
            ssFields = fieldnames(ss);
            for f=1:length(ssFields)
                ss.(ssFields{f}) = [ss.(ssFields{f}); nan(numRows,1)];
            end
        end
        
        function ss = cleanNaNs(ss)
            % Removes rows in which every entry is NaN (if only some of
            % them are NaN, it leaves the row intact)            
            ssFields = fieldnames(ss);
            L = dsSimpleStructUtil.length(ss);
            F = length(ssFields);
            nanFlag = false(L,F);
            
            % First loop through each field and flag the NaNs
            for f=1:F
                nanFlag(:,f) = isnan(ss.(ssFields{f}));
            end
            rowsToClean = all(nanFlag,2);
            
            % Now loop through each field and dlete the specified rows
            for f=1:F
                ss.(ssFields{f}) = ss.(ssFields{f})(~rowsToClean);
            end            
        end
        
        function ss = sortrows(ss,fieldToSortBy,howToSort)
            % Sort the struct by the values in fieldToSortBy.
            % howToSort = 'ascend' (default) or 'descend'
            if nargin < 3
                howToSort = 'ascend';
            end
            [~,pi_] = sort(ss.(fieldToSortBy), howToSort);
            
            ssFields = fieldnames(ss);
            for f=1:length(ssFields)
                ss.(ssFields{f}) = ss.(ssFields{f})(pi_);
            end
        end
        
    end        
    
end

