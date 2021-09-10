############################################################################
#
# This file provides a straightforward implementation of collection
# in pc groups, as a starting point for your own experiments.

import Base: inv # needed for `SyllableVector`

export abelianized,
       Collector,
       collector_from_pc_presentation,
       findfirst_uncollected_leftmost,
       freely_reduced,
       normalform,
       OscarInteger,
       Syllable,
       SyllableVector,
       test_collection


############################################################################
#
# define types, for convenience
#
# Syllable{T}   pair (generator, exponent) where exponent is of type `T`
# SyllableVector{T}

"""
    OscarInteger = Union{Base.Integer,fmpz}

In certain situations, we want to admit both Julia integers and `fmpz`.
"""
const OscarInteger = Union{Base.Integer,fmpz}

"""
    Syllable{T} = Tuple{Int,T}
    SyllableVector{T} = Vector{Syllable{T}}

The idea is to represent a word `x_{i1}^{e1}*x_{i2}^{e2}*...*x_{in}^{en}`,
where the `x_i` are group generators of type `Int`
and the `ei` are integers of type `T`,
by the sequence `[(i1, e1), (i2, e2), ..., (in, en)]`.
"""
const Syllable{T} = Tuple{Int,T}
const SyllableVector{T} = Vector{Syllable{T}}


############################################################################
#
# convert between `SyllableVector{T}` and `Vector{T}`

# `exponents` must contain the exponent of the `i`-th generator
# at position `i`
function SyllableVector{T}(exponents::Vector{T}) where T <: OscarInteger
    return Syllable{T}[(i, exponents[i]) for i in 1:length(exponents) if exponents[i] != 0]
end

# `syllables` must be in normal form
function Vector{T}(syllables::SyllableVector{T}, n::Int) where T <: OscarInteger
    v = zeros(T, n)
    for (i, e) in syllables
      v[i] = e
    end
    return v
end

# for convenience: compute the exponents w.r.t. a GAP pcgs
function SyllableVector{T}(pcgs::GAP.GapObj, pcelm::GAP.GapObj) where T <: OscarInteger
    exps = GAP.Globals.ExponentsOfPcElement(pcgs, pcelm)
    exps = Vector{T}(exps)
    return SyllableVector{T}(exps)
end

# for convenience: compute the syllable vector of an `FPGroupElem`
function SyllableVector{T}(fpelm::FPGroupElem) where T <: OscarInteger
    v = Vector{Int}(GAP.Globals.ExtRepOfObj(fpelm.X))
    return Syllable{T}[(v[i], v[i+1]) for i in 1:2:(length(v)-1)]
end

"""
    abelianized(v::SyllableVector{T}) where T <: OscarInteger

Return the normalized `SyllableVector{T}` obtained from abelianizing `v`,
that is, regarding the underlying generators as commuting.
"""
function abelianized(v::SyllableVector{T}) where T <: OscarInteger
    length(v) > 0 || return v
    return freely_reduced(sort(v))
end

############################################################################
#
# functions for `SyllableVector{T}`

"""
    freely_reduced(v::SyllableVector{T}) where T <: OscarInteger

Return a `SyllableVector{T}` that encodes the same group element as `v`
such that adjacent syllables belong to different generators.
"""
function freely_reduced(v::SyllableVector{T}) where T <: OscarInteger
    w = Syllable{T}[]
    for x in v
      if x[2] != 0
        if length(w) == 0
          push!(w, x)
        elseif w[end][1] == x[1]
          w[end] = (w[end][1], w[end][2] + x[2])
          if w[end][2] == 0
            pop!(w)
          end
        else
          push!(w, x)
        end
      end
    end
    return w
end

"""
    inv(v::SyllableVector{T}) where T <: OscarInteger

Return a `SyllableVector{T}` that encodes the inverse of `v`.
"""
function inv(v::SyllableVector{T}) where T <: OscarInteger
    return Syllable{T}[(v[i][1], -v[i][2]) for i in length(v):-1:1]
end


############################################################################
#
# define Collector{T}

"""
    Collector{T <: OscarInteger}

Objects of this type contain data needed to compute the normal forms of
words in pc groups that are given as `SyllableVector{T}`.

The fields are
- `findfirst_uncollected(v::SyllableVector{T})`,
  a function that returns an integer `i` such that either
  - the `i`-th syllable of `v` is an uncollected subword or
  - the `i`-th syllable and the first letter of the `(i+1)`-st syllable
    form an uncollected subword or
  - `0` if `v` is collected.

- `powers::Vector{T}`, where the `i`-th entry is the relative order `s_i`
  if `s_i` is finite, and `0` otherwise.

- `rules::Matrix{SyllableVector{T}}`, where the entry at `(i,j)` is
  `R_{i,j} = x_j^{-1} x_i x_j` if `1 <= j < i <= n`,
  `R_{i,j} = x_i x_j x_i^{-1}` if `1 <= i < j <= n`,
  `R_{i,i} = x_i^{s_i}` if `s_i` is finite.
"""
struct Collector{T <: OscarInteger}
    findfirst_uncollected::Function
    powers::Vector{T}
    rules::Matrix{SyllableVector{T}}
end


"""
    normalform(coll::Collector{T}, v::SyllableVector{T}; verbose::Bool = false) where T <: OscarInteger

Return a `SyllableVector{T}` that is the normal form of `v`,
w.r.t. collection given by the rules in `coll`.

If `verbose` is set to `true` then one line is printed per collection step,
showing the current word.
If the terminal capabilities admit then the (first) syllable of the
uncollected subword to be handled is highlighted in blue.
"""
function normalform(coll::Collector{T}, v::SyllableVector{T}; verbose::Bool = false) where T <: OscarInteger
    v = freely_reduced(v)
    i = coll.findfirst_uncollected(coll, v)
    while i > 0
      if verbose
        # Show the current word, try to highlight the `i-th syllable.
        println("# ", v[1:(i-1)]...,
                "\e[1m\e[38;2;0;0;255;249m", v[i], "\e[0m", v[(i+1):end]...)
      end
      w = v[1:i-1]
      exp = v[i][2]
      ii = v[i][1]
      s = coll.powers[ii]
      if s != 0 && (exp >= s || exp < 0)
        # Note that s cannot be negative.
        # write a = b*s + r, replace x_i^a by x_i^r * R_{i,i}^b
        # (Both b and r can be negative.)
        (b, r) = divrem(exp, s)
        if r < 0
          r = r + s
          b = b - 1
        end
        if r != 0
          push!(w, (v[i][1], r))
        end
        rule = coll.rules[ii, ii]
        if b < 0
          rule = inv(rule)
          b = -b
        end
        for j in 1:b
          append!(w, rule)
        end
        append!(w, v[(i+1):end])
      else
        # Which conjugator rule must be applied,
        # the one for x^a y or the one for x^a y^{-1}?
        x = v[i][1]
        a = v[i][2]
        y = v[i+1][1]
        b = v[i+1][2]
        @assert x > y "invalid minimal uncollected word"
        # We have either x_i^a x_j^b or x_i^a x_j^-b, with b > 0.
        if b > 0
          # rule below the diagonal
          rule = coll.rules[x, y]
          push!(w, (y, 1))
        else
          # rule above the diagonal
          rule = coll.rules[y, x]
          push!(w, (y, -1))
        end
        if a < 0
          rule = inv(rule)
          a = -a
        end
        for j in 1:a
          append!(w, rule)
        end
        if b > 1
          push!(w, (y, b-1))
        elseif b < -1
          push!(w, (y, b+1))
        end
        append!(w, v[(i+2):end])
      end
      v = freely_reduced(w)
      i = coll.findfirst_uncollected(coll, v)
    end
    return v
end


"""
    findfirst_uncollected_leftmost(coll::Collector{T}, v::SyllableVector{T}) where T <: OscarInteger

Return the leftmost position `i` such that either the `i`-th syllable of `v`
or this syllable together with the first letter of the `i+1`-st syllable
is a minimal uncollected subword.
If `v` is collected then `0` is returned.

`v` is assumed to be freely reduced, see [`freely_reduced`](@ref).
"""
function findfirst_uncollected_leftmost(coll::Collector{T}, v::SyllableVector{T}) where T <: OscarInteger
    for i in 1:length(v)
      s = coll.powers[v[i][1]]
      exp = v[i][2]
      if s != 0 && (exp >= s || exp < 0)
        return i
      end

      if i < length(v)
        if v[i][1] > v[i+1][1]
          return i
        end
      end
    end

    return 0
end


"""
    collector_from_pc_presentation(G::Oscar.GAPGroup, findfirst::Function = findfirst_uncollected_leftmost, ::Type{T} = Int) where T <: OscarInteger

Return a `Collector{T}` for the pc group `G`, such that `findfirst` is used
to find the first uncollected subword in a word.

(The data are just computed in GAP from the group `G.X`.)
"""
function collector_from_pc_presentation(G::Oscar.GAPGroup, findfirst::Function = findfirst_uncollected_leftmost, ::Type{T} = Int) where T <: OscarInteger
    @assert issolvable(G) "the group is not solvable"
    pcgs = GAP.Globals.Pcgs(G.X)
    n = length(pcgs)
    powers = Vector{T}(GAP.Globals.RelativeOrders(pcgs))
    rules = Base.Matrix{SyllableVector{T}}(undef, n, n)

    for i in 1:n
      # rules for powers
      rules[i, i] = SyllableVector{T}(pcgs, pcgs[i]^GAP.GapObj(powers[i]))

      # rules for conjugates
      for j in 1:(i-1)
        rules[i,j] = SyllableVector{T}(pcgs, pcgs[i]^pcgs[j])
        rules[j,i] = SyllableVector{T}(pcgs, pcgs[i]^(pcgs[j]^-1))
      end
    end

    return Collector{T}(findfirst, powers, rules)
end

# for convenience: enter a type but not a function
collector_from_pc_presentation(G::Oscar.GAPGroup, ::Type{T}) where T <: OscarInteger = collector_from_pc_presentation(G, findfirst_uncollected_leftmost, T)


"""
    test_collection(::Type{T}, G::Oscar.GAPGroup) where T <: OscarInteger

Run 100 tests for collection with a collector for `G` that is computed
by [`collector_from_pc_presentation`](@ref), comparing the collection result
for the concatenation of two random words with the product of the two
elements (computed by GAP in `G.X`),
and return `true` if the results coincide, and `false` otherwise.
"""
function test_collection(::Type{T}, G::Oscar.GAPGroup) where T <: OscarInteger
    coll = collector_from_pc_presentation(G, T)

    for i in 1:100
      g1 = GAP.Globals.Random(G.X)
      g2 = GAP.Globals.Random(G.X)
      prod = g1 * g2

      pcgs = GAP.Globals.Pcgs(G.X)
      exps1 = SyllableVector{T}(pcgs, g1)
      exps2 = SyllableVector{T}(pcgs, g2)
      exps = SyllableVector{T}(pcgs, prod)

      v = freely_reduced(vcat(exps1, exps2))
      nf = normalform(coll, v)
      if exps != nf
        println("$g1 * $g2 = $prod NOT $nf")
        return false
      end
    end
    return true
end
