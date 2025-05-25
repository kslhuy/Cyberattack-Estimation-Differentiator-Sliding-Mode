%% Dead zone function
function y = deadzone(x, dz)
    if abs(x) < dz
        y = 0;
    else
        y = x;
    end
end