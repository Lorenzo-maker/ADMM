function J = floatingBaseBodyJacZYX(q)

E = q(4:6);

C1 = cos(E(1));
C2 = cos(E(2));
C3 = cos(E(3));
S1 = sin(E(1));
S2 = sin(E(2));
S3 = sin(E(3));

J = [          C1*C2,          S1*C2,   -S2, 0, 0, 0;
     -S1*C3+C1*S2*S3, C1*C3+S1*S2*S3, C2*S3, 0, 0, 0;
      S1*S3+C1*S2*C3,-C1*S3+S1*S2*C3, C2*C3, 0, 0, 0;
      0, 0, 0,   -S2,  0, 1;
      0, 0, 0, C2*S3, C3, 0;
      0, 0, 0, C2*C3,-S3, 0];

end