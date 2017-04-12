function [value] = volume_step(x,a,thr)

if x>= thr
        value = a;
    else x< thr
        value = 0;
end


end


%%

