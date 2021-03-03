function y = pa2dbref(x)
    p0 = 1e-12;
    y = 10*log10(x/p0);
end