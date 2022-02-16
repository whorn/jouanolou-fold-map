Q = ZZ[x,y,z,x',y',z',u,v,u',v'] / ideal (x*(1-x)-y*z, x'*(1-x')-y'*z', x*u+y*v+x'*u'+y'*v'-1 ); -- JvJ model 

p = z*x'*y'*u*u'-x*x'*y'*v*u'-x*y'*z'*u'^2+z*y'^2*u*v'-x*y'^2*v*v'+2*x*x'*y'*u'*v'+x*y'^2*v'^2+x'*y'*v*u'+x*x'*u'^2+y'*z'*u'^2+y'^2*v*v'-2*x'*y'*u'*v'-y'^2*v'^2-z*y'*u+x*y'*v-x'*u'^2-y'*v+1;

ell = -z^2*y'*u^2-y*z*y'*v^2-2*z*x'*y'*v*u'-z*y'*z'*u'^2-2*z*y'^2*v*v'+2*z*x'*y'*u'*v'+z*y'^2*v'^2-2*z*y'*u*v+x*y'*v^2+z*x'*u*u'-x*x'*v*u'+z*x'*u'^2+z*y'*u*v'-x*y'*v*v'+2*z*y'*v-y'*v^2+x'*v*u'+y'*v*v'+z*u-x*v+v;

c = y'*(u*x+v*y)^2 + y*(u'*x'+v'*y')^2;

p' = -y*x'*z'*u'^2-x'*y'*z'*u'^2-2*y*y'*z'*u'*v'-2*y'^2*z'*u'*v'+y*x'*y'*v'^2+x'*y'^2*v'^2-y*y'*v'^2-y'^2*v'^2+2*y'*z'*u'-2*x'*y'*v'+2*y'*v'+x';

ell' = -y*z'^2*u'^2-y'*z'^2*u'^2+2*y*x'*z'*u'*v'+2*x'*y'*z'*u'*v'+y*y'*z'*v'^2+y'^2*z'*v'^2-2*y*z'*u'*v'-2*y'*z'*u'*v'+y*x'*v'^2+x'*y'*v'^2-2*x'*z'*u'-2*y'*z'*v'-y*v'^2-y'*v'^2+2*z'*u'-2*x'*v'+z'+2*v';

Fold = matrix{{p*p', c},{c *ell' * ell + p'^2 * ell + p^2 *ell', 1-p*p'}}; -- Fold map

R = ZZ[A,B,U,V,T]/ideal(A*U+V*B-1); -- ring needed for modeling degree 0 morphism 

deg0fold = sub(Fold, {x=>A*U, y =>B*U, z => A*V, x'=> (A-B)*(U+B), y' => -B*(U+B), z' => (A-B)*(A-B-U-V), u => 1+B*V, v =>-V^2, u' => 0, v' => -V^2 }); -- compute fold

M = matrix{{deg0fold_0_0 // (U*(U+B)), deg0fold_1_0//(U*(U+B))},{-deg0fold_0_1//(deg0fold_0_0 // (U*(U+B))),U*(U+B)}}; -- Turn J -> J map into a J -> A^2 - 0 map

s = (M_0_1+1)//M_1_1;

M'= M*matrix{{1,0},{-s,1}};

M'' = matrix{{1,0},{1,1}}*matrix{{1,(M'_0_0-1)},{0,1}}*M';

H = matrix{{1,0},{T*1,1}}*matrix{{1,T*(M'_0_0-1)},{0,1}}*M*matrix{{1,0},{-T*s,1}}; -- Final homotopy between fold and identity
