function X = hat(vec)

if size(vec,1)==3
    X = [      0,-vec(3), vec(2);
          vec(3),      0,-vec(1);
         -vec(2), vec(1),     0];
elseif size(vec,1)==6
    X = [hat(vec(4:6)),vec(1:3);zeros(4,1)];
end

end