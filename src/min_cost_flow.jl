#=
Primal-dual predictor corrector interior point solver for min-cost flow problems:

min c^Tx, subject to: B^Tx = b, x>=0, u>=x.

TODO: Addaptive termination.
TODO: Add adaptive regularization.
TODO: Active, non-active constraints
=#




function min_cost_flow{Tv,Ti}(mcfp::MCFproblem{Tv,Ti};
                              lapSolver = cholLap,
                              tol::Real=1e-6)

    edge_list = mcfp.edge_list
    m = size(edge_list,1)
    n = maximum(edge_list)
    B = sparse(collect(1:m), edge_list[:,1], 1.0, m, n) -
      sparse(collect(1:m), edge_list[:,2], 1.0, m, n)
    
    return min_cost_flow(B,
                         mcfp.costs,
                         mcfp.demands,
                         mcfp.capacities,
                         lapSolver = lapSolver,
                         tol = tol
                         )
end


function min_cost_flow{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                              c::Array{Tv,1},
                              b::Array{Tv,1},
                              u::Array{Tv,1};
                              #                              lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000)),
                              lapSolver = cholLap,
                              tol::Real=1e-6)
  # Problem dimensions.
  m = size(B)[1];
  n = size(B)[2];

  # Parameters for the algorithm.
  max_iter = 200;

  tol_p    = tol;
  tol_d    = tol;
  tol_gap  = tol;

  # Parameters for gap, residuals and Mehrotra's correction.
  eta   = 0;
  gamma = 1;

  # Compute initial point.
  (x,y,s,z) = ipm_min_cost_flow_initial_point(B,c,b,u,m,n,sddmSolver=cholSDDM);

  for k = 1:max_iter+1

    # Compute scaling parameters.
    w1 = sqrt(x./s);
    w2 = sqrt((u-x)./z);

    lambda1 = w1.*s;
    lambda2 = w2.*z;

    # Compute feasibility residuals and paramete mu.
    r_p = B'*x - b;
    r_d = c + B*y - s + z;
    mu  = (lambda1'*lambda1 + lambda2'*lambda2)/m;

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
    lambda1_sq = lambda1.*lambda1;
    lambda2_sq = lambda2.*lambda2;
    rhs_g1  = lambda1_sq;
    rhs_g2  = lambda2_sq;

      # note: w1, w2 and B do not change in here.
      # so, we can construct and reuse a solver
      
    (dx_a,dy_a,ds_a,dz_a) = ipm_directions_min_cost_flow(B,r_p,r_d,rhs_g1,rhs_g2,lambda1,lambda2,w1,w2,n,lapSolver=lapSolver);

    # Step-size and parameters.
    alpha_x = calstepsize(x,dx_a);
    alpha_x_up = calstepsize(u - x,-dx_a);
    alpha_s = calstepsize(s,ds_a);
    alpha_z = calstepsize(z,dz_a);
    alpha   = minimum([alpha_x;alpha_x_up;alpha_s;alpha_z]);

    rho   = ((x + alpha.*dx_a)'*(s + alpha.*ds_a) + (u - x - alpha.*dx_a)'*(z + alpha.*dz_a))/(x'*s + (u-x)'*z);
    sigma = (max(0,min(1,rho)))^3;

    # Combined direction.
    rhs_g1 = lambda1_sq - (sigma*mu).*ones(m) + gamma*(dx_a./w1).*(ds_a.*w1);
    rhs_g2 = lambda2_sq - (sigma*mu).*ones(m) - gamma*(dx_a./w2).*(dz_a.*w2);

    (dx,dy,ds,dz) = ipm_directions_min_cost_flow(B,(1-eta).*r_p,(1-eta).*r_d,
                                            rhs_g1,rhs_g2,lambda1,lambda2,w1,w2,n,lapSolver=lapSolver);

    # Update variables.
    alpha_x = calstepsize(x,dx);
    alpha_x_up = calstepsize(u - x,-dx);
    alpha_s = calstepsize(s,ds);
    alpha_z = calstepsize(z,dz);
    alpha   = minimum([alpha_x;alpha_x_up;alpha_s;alpha_z]);

    x = x + alpha.*dx;
    y = y + alpha.*dy;
    s = s + alpha.*ds;
    z = z + alpha.*dz;
  end
end

# Compute directions for nt_ipm.
function ipm_directions_min_cost_flow{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                                             rhs_p::Array{Tv,1},
                                             rhs_d::Array{Tv,1},
                                             rhs_g1::Array{Tv,1},
                                             rhs_g2::Array{Tv,1},
                                             lambda1::Array{Tv,1},
                                             lambda2::Array{Tv,1},
                                             w1::Array{Tv,1},
                                             w2::Array{Tv,1},
                                             m::Integer;
                                             lapSolver = cholLap
#                                             lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000))
    )

    Bt = B';
    
  # Solve saddle point system to compute directions dx,ds,dy.
  d1 = 1./(w1.*w1);
  d2 = 1./(w2.*w2);
  d  = 1./(d1 + d2);
  D = spdiagm(d);

  # Replacing this with the following faster code
  # L = B'*(D*B);
    # Adj = abs(spdiagm(diag(L)) - L)
    Adj = makeAdj(Bt,d)

  laInv = lapSolver((Adj+Adj')/2);

  lambda1_w1 = lambda1.*w1;
  lambda2_w2 = lambda2.*w2;

  rhs_comp = (-rhs_d - rhs_g1./lambda1_w1 + rhs_g2./lambda2_w2).*d;

  dy = laInv(rhs_p + Bt*rhs_comp);
  #dy = L\(rhs_p + Bt*rhs_comp);
  dx = -(B*dy).*d + rhs_comp;
  ds = -d1.*dx - rhs_g1./lambda1_w1;
  dz = d2.*dx - rhs_g2./lambda2_w2;

  return (dx,dy,ds,dz);
end

#=
Compute initial point by solving:

  maximize -b^Ty - u^Tz - 0.5y^Ty - 0.5*z^Tz
  subject to: c + By + z = 0

This is a regularized version of the dual of the min-cost flow problem by ignoring
the non-negativity constraints for the primal and the dual variables.

Then set s = x;

Then we keep the positive elements of the primal and dual variables x,y and we
further modify them so they are not too small or too large.
=#
function ipm_min_cost_flow_initial_point{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                                                c::Array{Tv,1},
                                                b::Array{Tv,1},
                                                u::Array{Tv,1},
                                                m::Integer,
                                                n::Integer;
                                                sddmSolver = cholSDDM
#                                                lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000))
)
  # Solve the optimization problem.
  L = B'*B + speye(n);

  laInv = sddmSolver(L);

  y = laInv(-b - B'*c + B'*u);
  #y = L\(-b - B'*c + B'*u);
  x = -c + B*y + u;
  z = x - u;
  s = x;

  # Keep only the positive elements.
  dx = max(-(3/2)*minimum(x),0);
  ds = max(-(3/2)*minimum(s),0);
  dz = max(-(3/2)*minimum(z),0);

  x = x + dx*ones(m);
  s = s + ds*ones(m);
  z = z + dz*ones(m);

  # Keep x <= u
  idx_large = x .>= u;
  if ~isempty(idx_large)
      x[idx_large] = (2/3)*u[idx_large];
  end

  # Modify them so they are not too small or too large.
  dx = 0.5*((x'*s + (u - x)'*z)/sum(s));
  ds = 0.5*((x'*s + (u - x)'*z)/sum(x));
  dz = 0.5*((x'*s + (u - x)'*z)/sum(u - x));

  x = x + dx.*ones(m);
  s = s + ds.*ones(m);
  z = z + dz.*ones(m);

  # Keep x <= u
  idx_large = x .>= u;
  if ~isempty(idx_large)
      x[idx_large] = (2/3)*u[idx_large];
  end

  return (x,y,s,z);
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

function makeAdj(Bt,w)
    n,m = size(Bt)
    bi = Bt.rowval[Bt.nzval.==1];
    bj = Bt.rowval[Bt.nzval.==-1];
    a = sparse([bj;bi],[bi;bj],[w;w],n,n)
    return a
end
