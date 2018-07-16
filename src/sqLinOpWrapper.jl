import Base.*
import Base.size
import Base.eltype
import LinearAlgebra.issymmetric

import LinearAlgebra.A_mul_B!

struct SqLinOp{Tv,Ti}
    issym::Bool
    value::Tv
    n::Ti
    multFn::Function #this is bad in terms of types?
   
    SqLinOp{Tv,Ti}(issym,value,n,multFn) where {Tv,Ti} = new(issym,value,n,multFn)
end

SqLinOp(issym,value::Tv,index::Ti,multFn) where {Tv,Ti} = SqLinOp{Tv,Ti}(issym,value,index,multFn)

eltype(A::SqLinOp{Tv,Ti}) where {Tv,Ti} = Tv

size(A::SqLinOp{Tv,Ti}, d::Ti) where {Tv,Ti} = A.n

size(A::SqLinOp{Tv,Ti}) where {Tv,Ti} = (A.n,A.n)

issymmetric(A::SqLinOp{Tv,Ti}) where {Tv,Ti} = A.issym

function *(A::SqLinOp{Tv,Ti}, b) where {Tv,Ti}
    return A.multFn(b)
end

function A_mul_B!(Y, A::SqLinOp{Tv,Ti}, B) where {Tv,Ti} 
    Y1 = A*B
    for i in 1:A.n
        Y[i] = Y1[i]
    end
end


function testId(n::Ti) where Ti
    return M = SqLinOp{Float64,Int64}(true,1.0,n,x -> x)
end

