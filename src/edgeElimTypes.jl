#=
  Types for the edgeElim solver
=#


"""
  LLp elements are all in the same column.
  row tells us the row, and val is the entry.
  val is set to zero for some edges that we should remove.
  next gives the next in the column.  It points to itself to terminate.
  reverse is the index into lles of the other copy of this edge,
  since every edge is stored twice as we do not know the order of elimination in advance.
"""
type LLp
    row::Int
    val::Float64
    next::LLp
    reverse::LLp

    LLp() = (x = new(0, 0.0); x.next = x; x.reverse = x)
    LLp(row, val, next, rev) = new(row, val, next, rev)
    LLp(row, val) = (x = new(row, val); x.next = x; x.reverse = x)
    LLp(row, val, next) = (x = new(row, val, next); x.reverse = x)
end

"""
  LLmatp is the data structure used to maintain the matrix during elimination.
  It stores the elements in each column in a singly linked list (only next ptrs)
  Each element is an LLp (linked list pointer).
  The head of each column is pointed to by cols.

  We probably can get rid of degs - as it is only used to store initial degrees.
"""
type LLmatp
    n::Int
    degs::Array{Int,1}
    cols::Array{LLp,1}
    lles::Array{LLp,1}
end

"""
  LLord: for an elimination of fixed order.
  Elements are all in the same column.
  Row tells us the row, and should be bigger than the column.
  Val is the entry.
  val is set to zero for some edges that we should remove.
  next gives the next in the column.  It points to itself to terminate.

"""
type LLord{Tind,Tval}
    row::Tind
    val::Tval
    next::LLord

    LLord{Tind,Tval}(row::Tind, val::Tval) = (x = new(row, val); x.next = x)
    LLord{Tind,Tval}(row::Tind, val::Tval, next) = (x = new(row, val, next))

end

"""
  LLMatOrd is the data structure used to maintain the matrix during elimination,
  when the order of elimination is pre-determined.
  It stores the elements in each column in a singly linked list (only next ptrs)
  Each element is an LLord (linked list pointer).
  The head of each column is pointed to by cols.

"""
type LLMatOrd{Tind,Tval}
    n::Tind
    cols::Array{LLord{Tind,Tval},1}
    lles::Array{LLord{Tind,Tval},1}
end





#=============================================================

LDLinv

=============================================================#

"""
  LDLinv contains the information needed to solve the Laplacian systems.
  It does it by applying Linv, then Dinv, then Linv (transpose).
  But, it is specially constructed for this particular solver.
  It does not explicitly make the matrix triangular.
  Rather, col[i] is the name of the ith col to be eliminated
"""
type LDLinv
    col::Array{Int,1}
    colptr::Array{Int,1}
    rowval::Array{Int,1}
    fval::Array{Float64,1}
    d::Array{Float64,1}
end

#=============================================================

EdgeElimPQ
the data strcture we use to keep track of degrees

=============================================================#

immutable EdgeElimPQElem
    prev::Int
    next::Int
    key::Int
end

"""
  An approximate priority queue.
  Items are bundled together into doubly-linked lists with all approximately the same key.
  minlist is the min list we know to be non-empty.
  It should always be a lower bound.
  keyMap maps keys to lists
"""
type EdgeElimPQ
    elems::Array{EdgeElimPQElem,1} # indexed by node name
    lists::Array{Int,1}
    minlist::Int
    nitems::Int
    n::Int
end
