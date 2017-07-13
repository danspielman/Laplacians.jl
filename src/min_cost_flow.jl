#=
Primal-dual predictor corrector interior point solver for min-cost flow problems:

min c^Tx, subject to: B^Tx = b, x>=0, u>=x.

TODO: Currently I am using backslash, check lines 131 - 136 and ipm_directions_min_cost_flow.
TODO: Some matvecs can be avoided by reusing matvecs from ipm_directions_min_cost_flow.
TODO: Comment out expensive printing statements, ctr + f to find them.

COMMENT: Centrality correctors have been commented out because they didn't seem to help.
COMMENT: I disable iterative refinement is only for dx and dy and not for dx, dy, ds and dz.
=#




function min_cost_flow{Tv,Ti}(mcfp::MCFproblem{Tv,Ti};
                              lapSolver = cholLap,
                              tol::Real=1e-6,
                              reg_p::Real=1e-10,
                              reg_d::Real=1e-10,
                              tol_ref = 1.0e-8)

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
                         tol = tol,
                         reg_p = reg_p,
                         reg_d = reg_d,
                         tol_ref = tol_ref)
end


function min_cost_flow{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                              c::Array{Tv,1},
                              b::Array{Tv,1},
                              u::Array{Tv,1};
                              lapSolver = cholLap,
                              tol::Real=1e-6,
                              reg_p::Real=1e-10,
                              reg_d::Real=1e-10,
                              tol_ref = 1.0e-8)
  # Problem dimensions.
  m = size(B)[1];
  n = size(B)[2];
    
  @printf("number of nodes =%d, number of edges=%d\n",n,m); 

  Bt = B'
    
  # Parameters for the algorithm.
  max_iter = 10000;

  tol_p    = tol;
  tol_d    = tol;
  tol_gap  = tol;
    
  # Parameters for gap, residuals and Mehrotra's correction.
  eta   = 0;
  gamma = 1;
  beta  = 1;
  beta2 = 1;
  alpha = 1;

  # Compute initial point.
  (x,y,s,z) = ipm_min_cost_flow_initial_point(B,c,b,u,m,n,sddmSolver=cholSDDM);
    
  # Dummy arrays
  ons = ones(m);
  zrs_m = zeros(m);
  zrs_n = zeros(n);
  x_old = x;
  y_old = y;

  for k = 1:max_iter+1

    # Compute feasibility residuals and paramete mu.
    r_p = Bt*x - b;
    r_d = c + B*y - s + z;
    mu  = (x'*s + (u-x)'*z)/(2*m);
    rel_gap = mu/(1 + abs(c'*x))
        
    max_x_s = maximum(x.*s);
    min_x_s = minimum(x.*s);
    @printf("maximum x.*s =%e, minimum x.*s=%e, mu=%e\n",max_x_s,min_x_s,mu[1]);  
        
    max_x_z = maximum((u-x).*z);
    min_x_z = minimum((u-x).*z);
    @printf("maximum (u-x).*z =%e, minimum (u-x).*z=%e, mu=%e\n",max_x_z,min_x_z,mu[1]);  

    # Check for termination.
    norm_r_p = norm(r_p)/(1 + norm(b));
    norm_r_d = norm(r_d)/(1 + norm(c));

    @printf("\nIteration %d, ||r_p||/||b||=%e, ||r_d||/||c||=%e, rel. gap=%e, alpha=%e\n", k, 
            norm_r_p, norm_r_d, rel_gap[1],alpha);

    if norm_r_p <= tol_p && norm_r_d <= tol_d && rel_gap[1] <= tol_gap
      println("Termination tolerance reached.");
      return (x,s,y);
    end

    if k > max_iter
      println("Maximum number of iteration reached.");
      return (x,s,y);
    end

    # Affine direction.
    rhs_g1 = x.*s;
    rhs_g2 = (u-x).*z;
        
    # Solve saddle point system to compute directions dx,ds,dy.        
    d1 = s./x;
    d2 = z./(u-x);
    theta = d1 + d2 + reg_p*ons;
    d  = 1./theta;
    D = spdiagm(d);     
        
    max_theta = maximum(theta);
    min_theta = minimum(theta);
    cond_theta = max_theta/min_theta;
    @printf("maximum theta =%e, minimum theta=%e, cond theta=%e\n", 
            max_theta,min_theta,cond_theta);
        
    max_theta_inv = maximum(d);
    min_theta_inv = minimum(d);
    cond_theta_inv = max_theta_inv/min_theta_inv;
    @printf("maximum theta inverse =%e, minimum theta inverse =%e, cond theta inverse=%e\n", 
            max_theta_inv,min_theta_inv,cond_theta_inv);

    # Compute adjacency and set lap solver.
    Adj = makeAdj(Bt,d);
    #laInv = lapSolver((Adj+Adj')/2);
    L = lap((Adj+Adj')/2);
    SDD = L + speye(n)*reg_d;
    laInv = SDD;
    
    # residuals for affine direction.
    r_p_tilde = r_p - reg_d*(y - y_old);
    r_d_tilde = r_d + reg_p*(x - x_old);
        
    rhs_d_saddle = -r_d_tilde - rhs_g1./x + rhs_g2./(u-x);
        
    # Compute affine direction.
    (dx_a,dy_a,ds_a,dz_a) = ipm_directions_min_cost_flow(B,Bt,d,d1,d2,theta,
                reg_p,reg_d,r_p_tilde,rhs_d_saddle,rhs_g1,rhs_g2,x,x_old,y,y_old,u,n,laInv);
        
    # Compute residuals for saddle point system of affine direction.
    res_p_saddle = Bt*dx_a - reg_d*dy_a + r_p_tilde;
    res_d_saddle = theta.*dx_a + B*dy_a - rhs_d_saddle;
        
    @printf("pred. norm(res_p_saddle_1) =%e, norm(res_d_saddle_1)=%e\n", 
            norm(res_p_saddle),norm(res_d_saddle));

    # Check if refinement is needed for affine direction and refine.
    if (norm(res_p_saddle) > tol_ref) || (norm(res_d_saddle) > tol_ref)          
            
        # Compute refined affine direction.
        rhs_normal = res_p_saddle + Bt*(res_d_saddle.*d);
        #dy = laInv(rhs_normal);
        dy_a_ref = laInv\rhs_normal;
        dx_a_ref = -(B*dy_a_ref).*d + res_d_saddle.*d;

        @printf("normal eq. residual =%e\n",norm(SDD*dy_a_ref - rhs_normal));  
            
        dx_a = dx_a + dx_a_ref;
        dy_a = dy_a + dy_a_ref;
            
        # Compute residuals for saddle point system of refined direction.
        res_p_saddle = Bt*dx_a - reg_d*dy_a + r_p_tilde;
        res_d_saddle = theta.*dx_a + B*dy_a - rhs_d_saddle;
            
        @printf("pred. norm(res_p_saddle_2) =%e, norm(res_d_saddle_2)=%e\n", 
            norm(res_p_saddle),norm(res_d_saddle));
            
        ds_a = -d1.*dx_a - rhs_g1./x;
        dz_a = d2.*dx_a - rhs_g2./(u-x);
    end
        
    alpha_x_a = calstepsize(x,dx_a);
    alpha_x_up_a = calstepsize(u - x,-dx_a);
    alpha_s_a = calstepsize(s,ds_a);
    alpha_z_a = calstepsize(z,dz_a);
    alpha_a = 0.99*minimum([alpha_x_a;alpha_x_up_a;alpha_s_a;alpha_z_a]);        
        
    x_a = x + alpha_a.*dx_a;
    y_a = y + alpha_a.*dy_a;
    s_a = s + alpha_a.*ds_a;
    z_a = z + alpha_a.*dz_a;  
        
    # Compute residual for affine solutions.
    r_p_a = Bt*x_a - b;
    r_d_a = c + B*y_a - s_a + z_a;   
        
    rho   = (x_a'*s_a + (u - x_a)'*z_a)/(x'*s + (u-x)'*z);
    sigma = (max(0,min(1,rho)))^3;
    mu_target = sigma*mu;
    #eta = 1 - sigma;
        
    x_a_s_a = x_a.*s_a;
    x_a_z_a = (u - x_a).*z_a;
    
    # Compute centrality corrector, i.e., correct just the outliers.    
    #rhs_g1 = zeros(m); 
    #rhs_g2 = zeros(m); 
     
    #idx_x_a_s_a_small = x_a_s_a .< (beta*mu_target[1]*ons);
    #idx_x_a_s_a_large = x_a_s_a .> ((mu_target[1]/beta)*ons);
    #idx_x_a_s_a_rest  = convert(Array{Bool},1 - (idx_x_a_s_a_small + idx_x_a_s_a_large));
        
    #idx_x_a_z_a_small = x_a_z_a .< (beta*mu_target[1]*ons);
    #idx_x_a_z_a_large = x_a_z_a .> ((mu_target[1]/beta)*ons);
    #idx_x_a_z_a_rest  = convert(Array{Bool},1 - (idx_x_a_z_a_small + idx_x_a_z_a_large));
        
    #rhs_g1[idx_x_a_s_a_small] = x_a_s_a[idx_x_a_s_a_small] - (mu_target[1]*beta);
    #rhs_g1[idx_x_a_s_a_large] = x_a_s_a[idx_x_a_s_a_large] - (mu_target[1]/beta);
    #rhs_g1[idx_x_a_z_a_rest]  = (x_a_s_a[idx_x_a_z_a_rest] - mu_target[1]).*beta2;
        
    #rhs_g2[idx_x_a_z_a_small] = x_a_z_a[idx_x_a_z_a_small] - (mu_target[1]*beta);
    #rhs_g2[idx_x_a_z_a_large] = x_a_z_a[idx_x_a_z_a_large] - (mu_target[1]/beta);
    #rhs_g2[idx_x_a_z_a_rest]  = (x_a_z_a[idx_x_a_z_a_rest] - mu_target[1]).*beta2;
           
    rhs_g1 = x.*(s_a - s) + s.*(x_a - x) + x.*s - mu_target.*ons + (x_a - x).*(s_a - s);
    rhs_g2 = (u-x).*(z_a - z) - z.*(x_a - x) + (u-x).*z - mu_target.*ons - (x_a - x).*(z_a - z);

    # residuals for saddle point for corrector direction.
    r_p_tilde = (1-eta[1]).*r_p_a - reg_d*(y - y_old);
    r_d_tilde = (1-eta[1]).*r_d_a + reg_p*(x - x_old);
        
    rhs_d_saddle = -r_d_tilde - rhs_g1./x + rhs_g2./(u-x);
        
    # Compute corrector direction.
    (dx,dy,ds,dz) = ipm_directions_min_cost_flow(B,Bt,d,d1,d2,theta,
                reg_p,reg_d,r_p_tilde,rhs_d_saddle,rhs_g1,rhs_g2,x,x_old,y,y_old,u,n,laInv);
        
    # Compute residuals for saddle point system of corrector direction.
    res_p_saddle = Bt*dx - reg_d*dy + r_p_tilde;
    res_d_saddle = theta.*dx + B*dy - rhs_d_saddle;
        
    @printf("corr. norm(res_p_saddle_1) =%e, norm(res_d_saddle_1)=%e\n", 
            norm(res_p_saddle),norm(res_d_saddle));

    # Check if refinement is needed for corrector direction and refine.
    if (norm(res_p_saddle) > tol_ref) || (norm(res_d_saddle) > tol_ref)

        # Compute refined corrector direction.
        rhs_normal = res_p_saddle + Bt*(res_d_saddle.*d);
        #dy = laInv(rhs_normal);
        dy_ref = laInv\rhs_normal;
        dx_ref = -(B*dy_ref).*d + res_d_saddle.*d;

        @printf("normal eq. residual =%e\n",norm(SDD*dy_ref - rhs_normal));  
            
        dx = dx + dx_ref;
        dy = dy + dy_ref;
            
        # Compute residuals for saddle point system of refined direction.
        res_p_saddle = Bt*dx - reg_d*dy + r_p_tilde;
        res_d_saddle = theta.*dx + B*dy - rhs_d_saddle;
            
        @printf("corr. norm(res_p_saddle_2) =%e, norm(res_d_saddle_2)=%e\n", 
            norm(res_p_saddle),norm(res_d_saddle));
            
        ds = -d1.*dx - rhs_g1./x;
        dz = d2.*dx - rhs_g2./(u-x);
    end       
        
    dx = dx + gamma.*dx_a  
    dy = dy + gamma.*dy_a  
    ds = ds + gamma.*ds_a  
    dz = dz + gamma.*dz_a
        
    # Update variables.
    alpha_x = calstepsize(x,dx);
    alpha_x_up = calstepsize(u - x,-dx);
    alpha_s = calstepsize(s,ds);
    alpha_z = calstepsize(z,dz);
    alpha = 0.99*minimum([alpha_x;alpha_x_up;alpha_s;alpha_z]);

    x_old = x;
    y_old = y;
        
    x = x + alpha.*dx;
    y = y + alpha.*dy;
    s = s + alpha.*ds;
    z = z + alpha.*dz;
  end
end

# Compute directions for nt_ipm.
function ipm_directions_min_cost_flow{Tv,Ti}(B::SparseMatrixCSC{Tv,Ti},
                                             Bt::SparseMatrixCSC{Tv,Ti},
                                             d, d1, d2, theta, reg_p, reg_d,
                                             rhs_p::Array{Tv,1},
                                             rhs_d_saddle::Array{Tv,1},
                                             rhs_g1::Array{Tv,1},
                                             rhs_g2::Array{Tv,1},
                                             x::Array{Tv,1},
                                             x_old::Array{Tv,1},
                                             y::Array{Tv,1},
                                             y_old::Array{Tv,1},
                                             u::Array{Tv,1},
                                             m::Integer,
                                             laInv
    )

  D = spdiagm(d);

  rhs_normal = rhs_p + Bt*(rhs_d_saddle.*d);
  #dy = laInv(rhs_normal);
  dy = laInv\rhs_normal;
  dx = -(B*dy).*d + rhs_d_saddle.*d;
  ds = -d1.*dx - rhs_g1./x;
  dz = d2.*dx - rhs_g2./(u-x);
    
  # XXX comment this out after testing.  
  @printf("normal eq. residual =%e\n",norm(B'*(D*(B*dy)) + reg_d*dy - rhs_normal));

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
