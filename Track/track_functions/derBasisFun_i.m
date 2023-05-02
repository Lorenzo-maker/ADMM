function ders = derBasisFun_i(u, i, p, U, derOrder)
ndu = zeros(p + 1, p + 1, class(u));
left = zeros(1, p + 1, class(u));  % initialize the left and right vectors (p. 69 NURBS book)
right = zeros(1, p + 1, class(u));
ndu(1,1) = 1;
for j = 1:p
   left(j) = u - U(i+1-j);
   right(j) = U(i+j) - u;
   saved = 0;
   for r = 0:j-1
       ndu(j+1, r+1) = right(r+1) + left(j-r);
       temp = ndu(r+1, j)./ndu(j+1,r+1);
       
       ndu(r+1, j+1) = saved + right(r+1).*temp;
       saved = left(j-r).*temp;
   end
   ndu(j+1,j+1) = saved;
end

ders = zeros(derOrder+1, p+1, class(u));
ders(1, :) = ndu(:, p+1);
a = zeros(2, p+1, class(u));
for r = 0 : p
   s1 = 0;
   s2 = 1;
   a(1,1) = 1;
   
   for k = 1:derOrder
       
       d = 0;
       rk = r - k;
       pk = p - k;
       
       if r >= k
          a(s2+1, 1) = a(s1+1, 1)./ndu(pk+2,rk+1);
          d = a(s2+1, 1).*ndu(rk+1, pk+1);
       end
       
       if rk >= -1; j1 = 1; else; j1 = -rk; end
       if r-1 <= pk; j2 = k-1; else; j2 = p-r;end

       for j = j1:j2
           a(s2+1, j+1) = (a(s1+1, j+1) - a(s1+1,j))./ndu(pk + 2, rk + j+1);
           d = d + a(s2+1, j+1).*ndu(rk + j+1, pk+1);
       end
       if r <= pk
           a(s2+1, k+1) = -a(s1+1, k)./ndu(pk + 2, r+1);
           d = d + a(s2 + 1, k + 1).*ndu(r+1, pk+1);
       end
       ders(k + 1, r + 1) = d;
       j = s1; s1 = s2; s2 = j;
   end
   
   
end

   r1 = p;
   
   for k = 1:derOrder
       ders(k+1, :) = ders(k+1, :)*r1;
       r1 = r1.*(p-k);
   end
% *** !!! errore pag 73 NURBS book ultimo loop fuori dal loop principale

% aggiungiamo zeri alle basi nulle
len = length(U);
numCtrlPts = len-p-1;
endLen = numCtrlPts -(i);

zers1 = zeros(derOrder+1, endLen);

zers2 = zeros(derOrder+1, i-p-1);

ders = [zers2, ders, zers1]';