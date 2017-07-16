import Base.*
import Base.size
import Base.eltype
import Base.issymmetric

import Base.LinAlg.A_mul_B!

struct SqLinOp{Tv,Ti}
    issym::Bool
    value::Tv
    n::Ti
    multFn::Function #this is bad in terms of types?
   
    SqLinOp{Tv,Ti}(issym,value,n,multFn) where {Tv,Ti} = new(issym,value,n,multFn)
end

SqLinOp{Tv,Ti}(issym,value::Tv,index::Ti,multFn) = SqLinOp{Tv,Ti}(issym,value,index,multFn)

eltype{Tv,Ti}(A::SqLinOp{Tv,Ti}) = Tv

size{Tv,Ti}(A::SqLinOp{Tv,Ti}, d::Ti) = A.n

size{Tv,Ti}(A::SqLinOp{Tv,Ti}) = (A.n,A.n)

issymmetric{Tv,Ti}(A::SqLinOp{Tv,Ti}) = A.issym

function *{Tv,Ti}(A::SqLinOp{Tv,Ti}, b)
    return A.multFn(b)
end

function A_mul_B!{Tv,Ti}(Y, A::SqLinOp{Tv,Ti}, B) 
    Y1 = A*B
    for i in 1:A.n
        Y[i] = Y1[i]
    end
end


function testId{Ti}(n::Ti)
    return M = SqLinOp{Float64,Int64}(true,1.0,n,x -> x)
end

