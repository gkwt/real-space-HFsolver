function [eval, evect, etot] = HFsolver(n, H_elec, guess, dx, eps)
% Input number of discrete steps, electronic Hamiltonian, initial guess, 
% convergence requirement, and the function will return the final
% eigenvalue and eigenvector for the Fock operator. Also total energy is
% calculated.

diff = 100.;
lastE = 0;
scfstep = 0;

[i,j,k] = ind2sub([n,n,n],[1:n^3]); % indices used to calculate distances

while abs(diff) > eps

    % Coulomb operator (should depend on the states)
    coulArray = zeros(1,n^3);
    weight = abs(guess).^2; % a "guess" state to calculate 
    for q = 1:n^3
        invdist = ((sqrt((i-i(q)).^2 + (j-j(q)).^2 + (k-k(q)).^2)).^(-1));
        invdist(isinf(invdist)) = 100.; % use large value
        invdist = weight'.*invdist;
        coulArray(q) = sum( invdist ) * dx.^3;
    end
    coulmat = diag(coulArray); % coulomb matrix
    
    % Exchange operator (should also depend on the states)
    weight = repmat(guess,[1,n^3]).*repmat(guess,[1,n^3])'; 
    exmat = zeros(n^3,n^3);
    for q = 1:n^3
        invdist = ((sqrt((i-i(q)).^2 + (j-j(q)).^2 + (k-k(q)).^2)).^(-1));
        invdist(isinf(invdist)) = 100.; % use large value
        exmat(q,:) = invdist * dx.^3;
    end
    exmat = weight.*exmat; % exchange matrix
    
    % create the fock operator
    F = H_elec + 2.*coulmat - exmat;
    [vect, E] = eig(F);
    % sort by eigenvalues
    [~,perm] = sort(diag(E));
    E = diag(E(perm,perm));
%     vect = vect(:,perm) ./ norm(vect(:,perm));
    
    % calculate total energy
%     H_tot = H_elec + coulmat - 0.5*exmat;
%     [~,E_hf] = eig(H_tot);
%     etot = E_hf(orbit,orbit);
%     E_hf = diag(vect' * H_tot * vect);
%     filled_E_hf = E_hf(E_hf<0);
%     etot = sum(filled_E_hf(orbit))

    % version 2 calculate total energy
    % since m goes to N/2, will end up only using first orbital (exchange
    % and coulomb the same)
    % interaction terms
    orb = vect(:,1) ./ norm(vect(:,1)); % the new lowest energy state
    intE = 0.;
    for q = 1:n^3
        for l = 1:n^3
            invdist = ((sqrt((i(l)-i(q)).^2 + (j(l)-j(q)).^2 + (k(l)-k(q)).^2)).^(-1));
            invdist(isinf(invdist)) = 100.; % arbitrarily large
            intE = intE + (orb(q)^2*orb(l)^2)*invdist;
        end
    end
    elecE = diag(ctranspose(orb)* H_elec*orb); % H_elec term
    
    etot = elecE + 0.5*intE

    % update parameters
    diff = abs( etot - lastE );
    guess = orb;
    lastE = etot;
    scfstep = scfstep+1;
    
end

eval = E; evect = guess;