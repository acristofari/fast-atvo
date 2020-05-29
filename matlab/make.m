function make()
    if (~verLessThan('matlab','9.4'))
        mex -R2018a -output fast_atvo fast_atvo_matlab.cpp ../fast_atvo.cpp
    else
        mex -largeArrayDims -output fast_atvo fast_atvo_matlab.cpp ../fast_atvo.cpp
    end
end