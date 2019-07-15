using JLD2

@show ARGS[1]
@show ARGS[2]
@show typeof(ARGS[2])

s = ARGS[2]

@show isa(s,AbstractString)

@load s dic
