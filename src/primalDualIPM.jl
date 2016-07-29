#=
Written by Anup Rao.
Primal-dual interior point method for a generic linear program of the form min_x c^Tx s.t. Ax = b, x>=0. This assumes the linear program is feasible and bounded and implements the practical implementation as outlined by Stephen Wright.

primalDualIPM takes as input a full row-rank matrix ‘A’, followed by vectors ‘c’ and 'b'. 
=#

function minCostFlow(Bt,b1,c1,u)

	#Bt = -dirEdgeVertexMat(a)
	#Bt = a
	m = size(Bt)[2]
	n = size(Bt)[1] 
	A = [Bt spzeros(n,m); speye(m) speye(m)]

	b = [b1;u]
	c = [c1;zeros(m,1)]

    #lapSolver=(H -> augTreeSolver(H,tol=1e-3,maxits=1000))
	lapSolver = (H -> lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000))

    #function lapSolver(H)
    #    return lapWrapSolver(augTreeSolver,H,tol=1e-8,maxits=1000)
    #end
    #f = lapSolver(H)
    #f(x1)
    #f(x2)



    shur = ((rhs,x,s) -> shurSolve(Bt,A,rhs,x,s,m,n,lapSolver))

	(x,y,s) = primalDualIPM(A,b,c,shur)

	return (x,y,s)



end



function primalDualIPM(A,b,c,shur)

tol = 1e-6
maxIter = 200
n1 = size(A)[1]
n2 = size(A)[2]
trunc = 0.995;

#solve = ((rb,rc,rd,x,s) -> solver1(A,rb,rc,rd,x,s,n1,n2))
solve = ((rb,rc,rd,x,s) -> solver3(Bt,A,rb,rc,rd,x,s,shur))

#computing the intial point
(x,y,s) = initialPoint1(A,b,c,shur)

for i=1:maxIter
   

mu_cur = BLAS.dot(x,s)/n2  #current mu
    
    
    
#predictor step
    rd = (x.*s)
    rb = (A*x - b)
    rc = A'*y + s - c

    if(norm([rb;rc;rd])<tol)
        display(i)
        break 
    end

	# d = (s.*(x.^-1))  #might need to add back the 1e+16 term
 #    D = diagm(d[:,1])

 #    H = [D -A';-A zeros(n1,n1)] 

 #    b1 = rc - (x.^-1).*rd
 #    b2 = rb



 #    dxy = H\[b1;b2] # need to substitute this with a solver
 #    dx = dxy[1:n2] 
 #    dy = dxy[n2+1:n1+n2]
 #    ds = -d.*dx - (x.^-1).*rd


    (dx,dy,ds) = solve(rb,rc,rd,x,s)
    
    id = find(dx.<0)
    if(isempty(id))
    	alpha_p=1
	else
		alpha_p = minimum(-x[id]./dx[id])
		alpha_p = min(1.0,alpha_p)
	end
	


	id = find(ds .< 0)
	if(isempty(id))
    	alpha_d=1
	else
		alpha_d = minimum(-s[id]./ds[id]);
		alpha_d = min(1,alpha_d);
	end

#compute corrector step parameters
	mu_aff = (x + alpha_p*dx)'*(s+alpha_d*ds)/n2
	sig = (mu_aff/mu_cur)^3

#corrector step
	rd = rd + dx.*ds - (sig*mu_cur)[1,1]
	(dx,dy,ds) = solve(rb,rc,rd,x,s)
    
    # b1 = rc - (x.^-1).*rd
    # b2 = rb



    # dxy = H\[b1;b2] # need to substitute this with a solver
    # dx = dxy[1:n2] 
    # dy = dxy[n2+1:n1+n2]
    # ds = -d.*dx - (x.^-1).*rd 

#compute the step size
	id = find(dx.<0)
	if(isempty(id))
    	alpha_pk=1
	else
		alpha_pk = minimum(-x[id]./dx[id]) 
		alpha_pk = min(1,trunc*alpha_pk)
	end
	

	id = find(ds .< 0)
	if(isempty(id)) 
    	alpha_dk=1
	else
		alpha_dk = minimum(-s[id]./ds[id])
		alpha_dk = min(1,trunc*alpha_dk)
	end

	x = x + alpha_pk*dx
	y = y + alpha_dk*dy
	s = s + alpha_dk*ds

end

return(x,y,s)

end

function initialPoint(A,b,c,shur)
A2 = A*A'
x0 = A'*(A2\b)
y0 = A2\(A*c)
s0 = c - A'*y0
del_x = max(-1.5*minimum(x0),0)
    del_s = max(-1.5*minimum(s0),0)
x0 = x0 + del_x
s0 = s0 + del_s

del_x = (0.5*BLAS.dot(x0,s0)/sum(s0))
del_s = 0.5*BLAS.dot(x0,s0)/sum(x0)


x0 = x0 + del_x
s0 = s0 + del_s

return (x0,y0,s0)

end

function initialPoint1(A,b,c,shur)
#A2 = A*A'
#x0 = A'*(A2\b)
n2 = size(A)[2]
x0 = A'*shur(b,ones(n2,1),ones(n2,1))
#y0 = A2\(A*c)
y0 = shur(A*c,ones(n2,1),ones(n2,1))
s0 = c - A'*y0
del_x = max(-1.5*minimum(x0),0)
    del_s = max(-1.5*minimum(s0),0)
x0 = x0 + del_x
s0 = s0 + del_s

del_x = (0.5*BLAS.dot(x0,s0)/sum(s0))
del_s = 0.5*BLAS.dot(x0,s0)/sum(x0)


x0 = x0 + del_x
s0 = s0 + del_s

return (x0,y0,s0)

end

function solver1(A,rb,rc,rd,x,s,n1,n2)
	d = (s.*(x.^-1))  #might need to add back the 1e+16 term
    D = spdiagm(d[:,1])

    H = [D -A';-A spzeros(n1,n1)] 

    b1 = rc - (x.^-1).*rd
    b2 = rb



    dxy = H\[b1;b2] # need to substitute this with a solver
    dx = dxy[1:n2] 
    dy = dxy[n2+1:n1+n2]
    ds = -d.*dx - (x.^-1).*rd

return(dx,dy,ds)

end

function solver2(A,rb,rc,rd,x,s,n1,n2)
	d = (s.^-1.*(x))  #might need to add back the 1e+16 term
    D = spdiagm(d[:,1])

    #H = [D -A';-A zeros(n1,n1)] 
    B = A*D*A'

    #b1 = rc - (x.^-1).*rd
    #b2 = rb

    Sinv = spdiagm((s.^-1)[:,1])
    X = spdiagm(x[:,1])

    dy = B\(-rb - A*(X*(Sinv*rc)) + A*(Sinv*rd)) 
    ds = -rc - A'*dy
    dx = -Sinv*rd - X*(Sinv*ds)
    #dxy = H\[b1;b2] # need to substitute this with a solver
    #dx = dxy[1:n2] 
    #dy = dxy[n2+1:n1+n2]
    #ds = -d.*dx - (x.^-1).*rd

return(dx,dy,ds)

end


function solver3(Bt,A,rb,rc,rd,x,s,shur)
    
    #x1 = x[1:m,1]
    #x2 = x[m+1:2*m,1]
    #s1 = s[1:m,1]
    #s2 = s[m+1:2*m,1]

    #d1 = (s1.^-1.*(x1))  #might need to add back the 1e+16 term
    #D1 = spdiagm(d1[:,1])

    #d2 = (s2.^-1.*(x2))  #might need to add back the 1e+16 term
    #D2 = spdiagm(d2[:,1])
    
    #wt = 1./(1./d1 + 1./d2)
    #la = Bt*spdiagm(wt[:,1])*Bt'        # this is the Shur complement
    
    #dinv = 1./(d1 + d2)
    #Dinv = spdiagm(dinv[:,1])

##
    Sinv = spdiagm((s.^-1)[:,1])
    X = spdiagm(x[:,1])
    rhs = -rb - A*(X*(Sinv*rc)) + A*(Sinv*rd) 
##

    rhs1 = rhs[1:n,1]
    rhs2 = rhs[n+1:end,1]

    #display(minimum(eigs(la)))
    #df = la\(rhs1 - Bt*D1*Dinv*rhs2)
    #laInv = lapSolver(la)
    #laInv = lapWrapSolver(augTreeSolver,la,tol=1e-6,maxits=1000)
    #df = laInv((rhs1 - Bt*D1*Dinv*rhs2))
    #dw = Dinv*(rhs2 - D1*Bt'*df)
##
   # dy = [ df;dw]
    #dy = shurSolve(Bt,A,rhs,x,s,m,n,lapSolver)
    dy = shur(rhs,x,s)
    ds = -rc - A'*dy
    dx = -Sinv*rd - X*(Sinv*ds) 
##

return(dx,dy,ds)
end

function shurSolve(Bt,A,rhs,x,s,m,n,lapSolver)

    x1 = x[1:m,1]
    x2 = x[m+1:2*m,1]
    s1 = s[1:m,1]
    s2 = s[m+1:2*m,1]

    d1 = (s1.^-1.*(x1))  #might need to add back the 1e+16 term
    D1 = spdiagm(d1[:,1])

    d2 = (s2.^-1.*(x2))  #might need to add back the 1e+16 term
    D2 = spdiagm(d2[:,1])
    
    wt = 1./(1./d1 + 1./d2)
    la = Bt*spdiagm(wt[:,1])*Bt'        # this is the Shur complement
    
    dinv = 1./(d1 + d2)
    Dinv = spdiagm(dinv[:,1])

##
    Sinv = spdiagm((s.^-1)[:,1])
    X = spdiagm(x[:,1])

    rhs1 = rhs[1:n,1]
    rhs2 = rhs[n+1:end,1]

    #display(minimum(eigs(la)))
    #df = la\(rhs1 - Bt*D1*Dinv*rhs2)
    laInv = lapSolver(la)
    df = laInv((rhs1 - Bt*D1*Dinv*rhs2))
    dw = Dinv*(rhs2 - D1*Bt'*df)
## 
    dy = [ df;dw]

return dy
end

