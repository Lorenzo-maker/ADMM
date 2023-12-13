function R = rotZYX(E)

C1 = cos(E(1));
C2 = cos(E(2));
C3 = cos(E(3));
S1 = sin(E(1));
S2 = sin(E(2));
S3 = sin(E(3));

R = [          C1*C2,-S1*C3+C1*S2*S3, S1*S3+C1*S2*C3;
               S1*C2, C1*C3+S1*S2*S3,-C1*S3+S1*S2*C3;
                 -S2,          C2*S3,          C2*C3];

end