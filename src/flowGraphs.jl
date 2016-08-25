function eyeFlowGraph(k)

    a = []
    if k<=1
        e = [1 3; 2 3; 1 4; 2 4]
        a = sparse(e[:,1],e[:,2],1,4,4)
        a = a + a'
        return a
    end
  

    a0 = eyeFlowGraph(k-1)
    n0 = size(a0)[1]

    if k>2
        bundle = sparse(diagm((ones(k*(k-2))),k))
        bundle = bundle + bundle'
    else
        bundle = sparse(1:2,1:2,zeros(2))
    end

    nb = size(bundle)[1]

    a1 = [a0 spzeros(n0,2*nb); spzeros(nb,n0) bundle spzeros(nb,nb); spzeros(nb,n0+nb) bundle]

    n1 = size(a1)[1]

    a = [spzeros(2,n1+2); spzeros(n1,2) a1]
    display(size(a1))

    sb1 = 2+n0
    sb2 = 2+n0+nb

    a[1,4] = 1
    a[4,1] = 1

    a[1,sb1+collect(1:k)] = 1
    a[sb1+collect(1:k),1] = 1
    a[3,sb1+nb-k+collect(1:k)] = 1
    a[sb1+nb-k+collect(1:k),3] = 1

    a[2,3] = 1
    a[3,2] = 1

    a[[4],sb2+collect(1:k)] = 1
    a[sb2+collect(1:k),[4]] = 1
    a[[2],sb2+nb-k+collect(1:k)] = 1
    a[sb2+nb-k+collect(1:k), [2]] = 1

    return sparse(a)

end


function eyeFlowGraph1(k)

    a = []
    if k<=1
        e = [1 3; 2 3; 1 4; 2 4]
        a = sparse(e[:,1],e[:,2],1,4,4)
        a = a + a'
        return a
    end
  

    a0 = eyeFlowGraph1(k-1)
    display(size(a0))
    n0 = size(a0)[1]

    
    bundle = sparse(diagm((ones(k*(k-1))),k))
    bundle = bundle + bundle'
    

    nb = size(bundle)[1]
    
    a1 = [a0 spzeros(n0,2*nb); spzeros(nb,n0) bundle spzeros(nb,nb); spzeros(nb,n0+nb) bundle]

    n1 = size(a1)[1]

    a = [spzeros(2,n1+2); spzeros(n1,2) a1]
    #display(size(a1))

    sb1 = 2+n0
    sb2 = 2+n0+nb

    a[1,4] = 1
    a[4,1] = 1

    a[1,sb1+collect(1:k)] = 1
    a[sb1+collect(1:k),1] = 1
    a[3,sb1+nb-k+collect(1:k)] = 1
    a[sb1+nb-k+collect(1:k),3] = 1

    a[2,3] = 1
    a[3,2] = 1

    a[[4],sb2+collect(1:k)] = 1
    a[sb2+collect(1:k),[4]] = 1
    a[[2],sb2+nb-k+collect(1:k)] = 1
    a[sb2+nb-k+collect(1:k), [2]] = 1

    return sparse(a)

end