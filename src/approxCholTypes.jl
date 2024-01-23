#=
  Structs for the approxChol solver
=#


"""
  LLp elements are all in the same column.
  row tells us the row, and val is the entry.
  val is set to zero for some edges that we should remove.
  next gives the next in the column.  It points to itself to terminate.
  reverse is the index into lles of the other copy of this edge,
  since every edge is stored twice as we do not know the order of elimination in advance.
"""
mutable struct LLp{Tind,Tval}
    row::Tind
    val::Tval
    next::LLp{Tind,Tval}
    reverse::LLp{Tind,Tval}

    LLp{Tind,Tval}() where {Tind,Tval} = (x = new(zero(Tind), zero(Tval)); x.next = x; x.reverse = x)
    LLp{Tind,Tval}(row, val, next, rev) where {Tind,Tval} = new(row, val, next, rev)
    LLp{Tind,Tval}(row, val) where {Tind,Tval} = (x = new(row, val); x.next = x; x.reverse = x)
    LLp{Tind,Tval}(row, val, next) where {Tind,Tval} = (x = new(row, val, next); x.reverse = x)
end

"""
  LLmatp is the data structure used to maintain the matrix during elimination.
  It stores the elements in each column in a singly linked list (only next ptrs)
  Each element is an LLp (linked list pointer).
  The head of each column is pointed to by cols.

  We probably can get rid of degs - as it is only used to store initial degrees.
"""
mutable struct LLmatp{Tind,Tval}
    n::Int64
    degs::Array{Tind,1}
    cols::Array{LLp{Tind,Tval},1}
    lles::Array{LLp{Tind,Tval},1}
end

# these are the types we use with a fixed ordering
struct LLord{Tind,Tval}
    row::Tind
    next::Tind
    val::Tval
end

mutable struct LLMatOrd{Tind,Tval}
    n::Int64
    cols::Array{Tind,1}
    lles::Array{LLord{Tind,Tval},1}
end

struct LLcol{Tind,Tval}
      row::Tind
      ptr::Tind
      val::Tval
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
mutable struct LDLinv{Tind,Tval}
    col::Array{Tind,1}
    colptr::Array{Tind,1}
    rowval::Array{Tind,1}
    fval::Array{Tval,1}
    d::Array{Tval,1}
end

#=============================================================

ApproxCholPQ
the data strcture we use to keep track of degrees

=============================================================#

struct ApproxCholPQElem{Tind}
    prev::Tind
    next::Tind
    key::Tind
end

"""
  An approximate priority queue.
  Items are bundled together into doubly-linked lists with all approximately the same key.
  minlist is the min list we know to be non-empty.
  It should always be a lower bound.
  keyMap maps keys to lists
"""
mutable struct ApproxCholPQ{Tind}
    elems::Array{ApproxCholPQElem{Tind},1} # indexed by node name
    lists::Array{Tind,1}
    minlist::Int64
    nitems::Int64
    n::Int64
    # RG21 this is just a quick hack to make it work with multiedge
    split::Tind
    last_vertex::Union{Nothing,Tind}
end
