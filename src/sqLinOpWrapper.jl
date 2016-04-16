import Base.*
import Base.size
import Base.eltype

immutable SqLinOp{Tv,Ti}
    issym::Bool
    value::Tv
    n::Ti
    multFn::Function #this is bad in terms of types?
   
    SqLinOp(issym,value,n,multFn) = new(issym,value,n,multFn)
end

SqLinOp{Tv,Ti}(issym,value::Tv,index::Ti,multFn) = SqLinOp{Tv,Ti}(issym,value,index,multFn)

eltype{Tv,Ti}(A::SqLinOp{Tv,Ti}) = Tv

size{Tv,Ti}(A::SqLinOp{Tv,Ti}, d::Ti) = A.n

function *{Tv,Ti}(A::SqLinOp{Tv,Ti}, b::Array{Tv,1})
    return A.multFn(b)
end

function testId()
    return M = SqLinOp{Float64,Int64}(true,1.0,1,x -> x)
end
