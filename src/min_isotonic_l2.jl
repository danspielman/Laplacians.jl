#=
Primal-dual predictor corrector interior point solver for isotonic l2 problems:

min ||x-b||_2^2 subject to: Bx<=0.

TODO: Addaptive termination.
TODO: Add adaptive regularization.
TODO: Active, non-active constraints
=#

function min_isotonic_l2{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                                b::Array{Tv,1},
                                lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000)),
                                tol::Real=1e-6)
  # Problem dimensions.
  m = size(B)[1];
	n = size(B)[2];

  # Parameters for the algorithm.
  max_iter = 100;

  tol_p    = tol;
  tol_d    = tol;
  tol_gap  = tol;

  # Parameters for gap, residuals and Mehrotra's correction.
  eta   = 0;
  gamma = 1;

  # Compute initial point.
  (x,s,y) = ipm_l2_isotonic_initial_point(B,b,m,n,lapSolver);

  for k = 1:max_iter+1

    # Compute scaling parameters.
    w = sqrt(s./y);

    lambda = w.*y;

    # Compute feasibility residuals and paramete mu.
    r_p = B*x + s;
    r_d = x - b + B'*y;
    mu  = (lambda'*lambda)/m;

    # Check for termination.
    norm_r_p = norm(r_p);
    norm_r_d = norm(r_d);

    @printf("Iteration %d, ||r_p||=%f, ||r_d||=%f, mu=%f\n", k, norm_r_p, norm_r_d, mu[1]);

    if norm_r_p <= tol_p && norm_r_d <= tol_d && mu[1] <= tol_gap
      println("Termination tolerance reached.");
      return (x,s,y);
    end

    if k > max_iter
      println("Maximum number of iteration reached.");
      return (x,s,y);
    end

    # Affine direction.
    lambda_sq = lambda.*lambda;
    rhs_path  = lambda_sq;

    (dx_a,ds_a,dy_a) = ipm_directions_isotonic_l2(B,r_p,r_d,rhs_path,lambda,w,n,lapSolver);

    # Step-size and parameters.
    alpha_s = calstepsize(s,ds_a);
    alpha_y = calstepsize(y,dy_a);
    alpha   = min(alpha_s,alpha_y);

    rho   = ((s + alpha.*ds_a)'*(y + alpha.*dy_a))/(s'*y);
    sigma = (max(0,min(1,rho)))^3;

    # Combined direction.
    rhs_path = lambda_sq - (sigma*mu).*ones(m) + gamma*(ds_a./w).*(dy_a.*w);

    (dx,ds,dy) = ipm_directions_isotonic_l2(B,(1-eta).*r_p,(1-eta).*r_d,
                                            rhs_path,lambda,w,n,lapSolver);

    # Update variables.
    alpha_s = calstepsize(s,ds);
    alpha_y = calstepsize(y,dy);
    alpha   = min(alpha_s,alpha_y);

    x = x + alpha.*dx;
    s = s + alpha.*ds;
    y = y + alpha.*dy;
  end
end

# Compute directions for nt_ipm.
function ipm_directions_isotonic_l2{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                                           rhs_p::Array{Tv,1},
                                           rhs_d::Array{Tv,1},
                                           rhs_path::Array{Tv,1},
                                           lambda::Array{Tv,1},
                                           w::Array{Tv,1},
                                           n::Integer,
                                           lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000)))
  # Solve saddle point system to compute directions dx,ds,dy.
  d = 1./(w.*w);
  D = spdiagm(d);

  L = B'*(D*B) + speye(n);

  #laInv = lapSolver(L);

  lambda_w = lambda.*w;

  #dx = laInv(-rhs_d + B'*(rhs_path./lambda_w) - B'*(d.*rhs_p));
  dx = L\(-rhs_d + B'*(rhs_path./lambda_w) - B'*(d.*rhs_p));
  ds = -rhs_p - B*dx;
  dy = -d.*ds - rhs_path./lambda_w;

  return (dx,ds,dy);
end

#=
Compute initial point by solving:

  maximize -1/2 x'*x + -1/2 y'*y + b'*b
  subject to: x - b + B'y = 0

This is a regularized version of the dual of the l2-isotonic problem by ignoring
the non-negativity constraints for the primal and the dual variables.

Then set s = -Bx;

Then we keep the positive elements of the primal and dual variables x,y and we
further modify them so they are not too small or too large.
=#
function ipm_l2_isotonic_initial_point{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                                              b::Array{Tv,1},
                                              m::Integer,
                                              n::Integer,
                                              lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000)))
  # Solve the optimization problem.
  L = B'*B + speye(n);

  #laInv = lapSolver(L);

  #x = laInv(b);
  x = L\b;
  y = B*x;
  s = -y;

  # Keep only the positive elements.
  ds = max(-(3/2)*minimum(s),0);
  dy = max(-(3/2)*minimum(y),0);

  s = s + ds*ones(m);
  y = y + dy*ones(m);

  # Modify them so they are not too small or too large.
  ds = 0.5*((s'*y)/sum(y));
  dy = 0.5*((s'*y)/sum(s));

  s = s + ds.*ones(m);
  y = y + dy.*ones(m);

  return (x,s,y);
end

# Calculate step-size for positive orthant.
function calstepsize{Tv}(x::Array{Tv,1},dx::Array{Tv,1},maxstepsize::Float64 = 0.99)
  stepsizes = -x./dx;

  idx_pos = find(stepsizes .> 0);
  if isempty(idx_pos)
      return maxstepsize;
  else
      return minimum([0.999*stepsizes[idx_pos]; maxstepsize]);
  end
end
