function X = adjointStar(Gst)
%
% adjoint(Gst)
% trasformata aggiunta di Gst. Matrice 6x6 che può trasformare
% twist congruenti in base alla trasformazione omogenea di Gst.


R = Gst(1:3,1:3);
d = Gst(1:3,4);
X = [R,zeros(3,3);hat(d)*R,R];

end