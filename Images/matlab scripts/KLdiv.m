function KLdiv = KLdiv(p,q)

KLdiv = sum(p.*log(p./q));

end