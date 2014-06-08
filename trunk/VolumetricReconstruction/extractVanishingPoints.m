function [A_r b_r]=extractVanishingPoints(Hkk,Sub_cc)
Hkk = Sub_cc * Hkk;

% Extract vanishing points (direct and diagonals):

V_hori_pix = Hkk(:,1);
V_vert_pix = Hkk(:,2);
V_diag1_pix = (Hkk(:,1)+Hkk(:,2))/2;
V_diag2_pix = (Hkk(:,1)-Hkk(:,2))/2;

V_hori_pix = V_hori_pix/norm(V_hori_pix);
V_vert_pix = V_vert_pix/norm(V_vert_pix);
V_diag1_pix = V_diag1_pix/norm(V_diag1_pix);
V_diag2_pix = V_diag2_pix/norm(V_diag2_pix);

a1 = V_hori_pix(1);
b1 = V_hori_pix(2);
c1 = V_hori_pix(3);

a2 = V_vert_pix(1);
b2 = V_vert_pix(2);
c2 = V_vert_pix(3);

a3 = V_diag1_pix(1);
b3 = V_diag1_pix(2);
c3 = V_diag1_pix(3);

a4 = V_diag2_pix(1);
b4 = V_diag2_pix(2);
c4 = V_diag2_pix(3);

A_r = [a1*a2  b1*b2;
    a3*a4  b3*b4];

b_r = -[c1*c2;c3*c4];
end
