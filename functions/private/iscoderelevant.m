function useCode = iscoderelevant(str)

    useCode = true;
    if contains(str,'plot')
        useCode = false;
    end
    if contains(str,'clf')
        useCode = false;
    end
    if contains(str,'clear')
        useCode = false;
    end
    if contains(str,'validate')
        useCode = false;
    end
    if contains(str,'Parameters')
        useCode = false;
    end
    if contains(str,'function')
        useCode = false;
    end
    if contains(str,'clc')
        useCode = false;
    end
    if contains(str,'axis')
        useCode = false;
    end        
    
end