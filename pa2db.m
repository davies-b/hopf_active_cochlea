function y = pa2db(x)
    p0 = 2e-5;
    y = 20*log10(x/p0);
end