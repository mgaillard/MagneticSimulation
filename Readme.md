# Magnetic simulation
Axis symmetric simulation of magnetic field using finite differences.

## TODO
- Contour lines of the magnetic field
- Add screenshots to the Readme
- Output the plot of the magnetic field on the line r=0.0
- Add default visualizations in Paraview
- Split as a lib (Eigen-only), console App (Qt), and test project (Catch2).
- Add a test with a spire and compare to the theory formula
- Animation with interpolation of a parameter between two values
- Export the scene in 3D with an OBJ mesh
- Translate the code to english
- Improve and refactor code with C++11
- Change the smoothing convolution with Eigen convolution


## Testing code in Matlab
```matlab
clear

mu0 = 4.0 * pi * 1e-7;
Nfilt=20;

N = 600;
M = 600;
h = 1e-4;
k = 1e-4;
r = (0:h:(N-1)*h);
z = (-(M-2)/2*k:k:M/2*k);
r0 = 1e-2;

% Definition of permeability
mu = mu0 * ones(M, N);
eta0 = 1.0 ./ mu;

% [eta, dedr, dedz] = LissageSpatial(eta0, Nfilt, M, N, k, h);
% Fake lissage because mu is constant in our example
eta = eta0;
dedr = zeros(M, N);
dedz = zeros(M, N);

% Definition of density of currents
I = 1.0;
jd0 = zeros(M, N);
jd0(round(M / 2), r0 / h) = 1.0;
jd = jd0 * I / (h * k * sum(sum(jd0)));

% Convert matrices to vectors
Eta = Matr2Vect(eta, M, N);
Dedr = Matr2Vect(dedr, M, N);
Dedz = Matr2Vect(dedz, M, N);
J = Matr2Vect(jd, M, N);

D = CalculMatrice(N, M, Dedr, Dedz, h, k, Eta);

Al = Solution(Eta, D, J);
A = Vect2Matr(Al, M, N);

Bz = zeros(M,N);
Br = zeros(M,N);
for i=2:N-1
    for j=2:M-1
        Bz(j,i) = A(j,i)/r(i)+(A(j,i+1)-A(j,i-1))/(2*h);
        Br(j,i) = (-A(j+1,i)+A(j-1,i))/(2*k);
    end
end
Bz(:,1) = Bz(:,2);
Bt = sqrt(Br.^2 + Bz.^2);

% Check result versus theory
Bcentre = Bz(round(M/2), 1);
BcentreTh = mu0*I/(2*r0);
disp(Bcentre)
disp(BcentreTh)


Bth = (mu0 * I * r0 * r0)./(2.0 * power((z .* z) + (r0 * r0), 1.5));

plot(z, Bt(:, 1), z, Bth);
```