function rIonCoul = rIonCoul(n, Z, centre)
% Generate discrete matrix of 3D Coulomb interaction with ion. Position of
% the ion is at centre, with Z protons. 3D space is discretised into n
% parts. Gives restricted HF interaction, so there are two electrons.

centre = double(centre);
ionCoul = zeros(n,n,n);

for k = 1:n % cycle through 'z' coordinate
   for j = 1:n % cycle through 'y' coordinate
       for i = 1:n % cycle through 'x' coordinate
           ionCoul(i,j,k) = pdist([i,j,k;centre],'euclidean');
       end
   end
end

ionCoul = reshape(ionCoul,[1,n^3]);

ionCoul(find(~ionCoul)) = 0.001; % arbitrarily small number for 0s
ionCoul = -2.*Z*(ionCoul.^(-1)); % multiply by 2 for electrons, and Z for ions

rIonCoul = diag(ionCoul);
   





