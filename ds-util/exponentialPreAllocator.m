classdef exponentialPreAllocator < handle
    % A builder for datasets, structs, and simpleStructs to increase speed

    properties (Access = private)
        data;
        length;
    end

    methods
        function exPA = exponentialPreAllocator
            exPA.data = [];
            exPA.length = 0;
        end

        function addArrayOfDs(exPA, arrayOfDs)
            arrayOfSs = cell(size(arrayOfDs,1),1);
            for k=1:size(arrayOfDs,1),
                if ~isempty(arrayOfDs{k})
                    arrayOfSs{k} = dsSimpleStructUtil.dataset2SimpleStruct(arrayOfDs{k});
                end
            end
            exPA.addArrayOfSimpleStruct(arrayOfSs);
        end

        function addArrayOfStruct(exPA, arrayOfSt)
            arrayOfSs = cell(size(arrayOfSt,1),1);
            for k=1:size(arrayOfSt,1),
                if ~isempty(arrayOfSt{k})
                    arrayOfSs{k} = dsStructUtil.struct2SimpleStruct(arrayOfSt{k});
                end
            end
            exPA.addArrayOfSimpleStruct(arrayOfSs);
        end

        function addArrayOfSimpleStruct(exPA, arrayOfSs)
            d = exPA.data;
            len = exPA.length;
            
            for k=1:size(arrayOfSs,1),
                ss = arrayOfSs{k};
                if isempty(ss),
                    continue;
                end

                addLength = dsSimpleStructUtil.length(ss);

                if isempty(d),
                    exPA.addToEmpty(ss);
                    d = exPA.data;
                    len = exPA.length;
                    continue;
                end

                neededSize = len + addLength;
                if neededSize > dsSimpleStructUtil.length(d),
                    exPA.data = d;
                    exPA.increaseSize(neededSize);
                    d = exPA.data;
                end

                d = dsSimpleStructUtil.setSimpleStructSubset(d, len + (1:addLength), ss);
                len = len + addLength;
            end

            exPA.data = d;
            exPA.length = len;
        end
        function addDs(exPA, ds)
            if isempty(ds),
                return;
            end
            ss = dsSimpleStructUtil.dataset2SimpleStruct(ds);
            exPA.addSimpleStruct(ss);
        end

        function addStruct(exPA, st)
            if isempty(st),
                return;
            end
            ss = dsStructUtil.struct2SimpleStruct(st);
            exPA.addSimpleStruct(ss);
        end

        function addSimpleStruct(exPA, ss)
            if isempty(ss),
                return;
            end

            d = exPA.data;
            len = exPA.length;

            addLength = dsSimpleStructUtil.length(ss);

            if isempty(d),
                exPA.addToEmpty(ss);
                return;
            end

            neededSize = len + addLength;
            if neededSize > dsSimpleStructUtil.length(d),
                exPA.data = d;
                exPA.increaseSize(neededSize);
                d = exPA.data;
            end

            d = dsSimpleStructUtil.setSimpleStructSubset(d, len + (1:addLength), ss);
            len = len + addLength;

            exPA.data = d;
            exPA.length = len;
        end

        function ds = buildDs(exPA)
            if isempty(exPA.data)
                ds = [];
                return;
            end
            ss = exPA.buildSimpleStruct;
            ds = dsSimpleStructUtil.simpleStruct2Dataset(ss);
        end

        function st = buildStruct(exPA)
            if isempty(exPA.data)
                st = [];
                return;
            end
            ss = exPA.buildSimpleStruct;
            st = dsStructUtil.simpleStruct2Struct(ss);
        end

        function ss = buildSimpleStruct(exPA)
            if isempty(exPA.data)
                ss = [];
                return;
            end
            ss = dsSimpleStructUtil.getSimpleStructSubset(exPA.data, 1:exPA.length);
        end

    end

    methods (Access = private)

        function increaseSize(exPA, neededSize)
            while dsSimpleStructUtil.length(exPA.data) < neededSize,
                exPA.data = dsSimpleStructUtil.vertcat(exPA.data, exPA.data);
            end
        end

        function addToEmpty(exPA, ss)
            exPA.data = ss;
            exPA.length = dsSimpleStructUtil.length(ss);
        end

    end

end

