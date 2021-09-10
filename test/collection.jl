@testset "entered by hand: dihedral group of order 8" begin
  T = Int
  id = SyllableVector{T}()
  rules = Base.Matrix{SyllableVector{T}}(undef, 3, 3)

  powers = [2, 2, 2]

  rules[1,1] = [(2, 1)]
  rules[2,2] = id
  rules[3,3] = id

  rules[1,2] = [(2, 1)]
  rules[1,3] = [(2, 1), (3, 1)]
  rules[2,1] = [(2, 1)]
  rules[2,3] = [(3, 1)]
  rules[3,1] = [(2, 1), (3, 1)]
  rules[3,2] = [(3, 1)]
  
  coll = Collector{T}(findfirst_uncollected_leftmost, powers, rules)
  
  @test normalform(coll, [(1, 2)]) == [(2, 1)]
  @test normalform(coll, [(1, 4)]) == id
  @test normalform(coll, [(1, 1), (2, 1), (1, 1)]) == id
  @test normalform(coll, [(3, 1), (1, 1)]) == [(1, 1), (2, 1), (3, 1)]
end

@testset "entered by hand: dihedral group of order 16" begin
  # The presentation is shown on page 35 of Eamonn's talks:
  # < x1, x2, x3, x4 | x1^2 = 1, x2^2 = x3 x4,
  #                    x3^2 = x4, x4^2 = 1,
  #                    x1^-1 x2 x1 = x2 x3, x1^-1 x3 x1 = x3 x4,
  #                    x2^-1 x3 x2 = x3, x1^-1 x4 x1 = x4,
  #                    x2^-1 x4 x2 = x4, x3^-1 x4 x3 = x4 >
  T = Int8
  id = SyllableVector{T}()
  rules = Base.Matrix{SyllableVector{T}}(undef, 4, 4)

  powers = [2, 2, 2, 2]

  rules[1,1] = id
  rules[2,2] = [(3, 1), (4, 1)]
  rules[3,3] = [(4, 1)]
  rules[4,4] = id

  rules[1,2] = [(2, 1), (3, 1)]
  rules[1,3] = [(3, 1), (4, 1)]
  rules[1,4] = [(4, 1)]
  rules[2,1] = [(2, 1), (3, 1)]
  rules[2,3] = [(3, 1), (4, 1)]
  rules[2,4] = [(4, 1)]
  rules[3,1] = [(3, 1), (4, 1)]
  rules[3,2] = [(3, 1)]
  rules[3,4] = [(3, 1)]
  rules[4,1] = [(4, 1)]
  rules[4,2] = [(4, 1)]
  rules[4,3] = [(4, 1)]

  coll = Collector{T}(findfirst_uncollected_leftmost, powers, rules)

  @test normalform(coll, Tuple{Int64, T}[(3, 1), (2, 1), (1, 1)]) == [(1, 1), (2, 1)]
end

@testset "entered by hand: infinite dihedral group" begin
  T = fmpz
  id = SyllableVector{T}()

  # Let a = [-1 0; 0 1], b = [-1 -1; 0 1] (both of order 2);
  # a and b generate the infinite dihedral group.
  # Set x = a*b, it has infinite order.

  # 1. Take the polycyclic sequence [a, x].
  #    We have a^-1*x*a  = a*x*a^-1 = x^-1
  rules1 = Base.Matrix{SyllableVector{T}}(undef, 2, 2)

  powers1 = [2, 0]

  rules1[1,1] = id

  rules1[1,2] = [(2, -1)]
  rules1[2,1] = [(2, -1)]

  coll1 = Collector{Int}(findfirst_uncollected_leftmost, powers1, rules1)

  @test normalform(coll1, [(1, 1), (2, 1), (1, 1)]) == [(2, -1)]
  @test normalform(coll1, [(2, -1), (1, 1)]) == [(1, 1), (2, 1)]
  @test normalform(coll1, [(2, -2), (1, 1)]) == [(1, 1), (2, 2)]

  # 2. Take the polycyclic sequence [x, a].
  #    We have x^-1*a*x  = x^-2*a and x*a*x^-1 = x^2*a
  rules2 = Base.Matrix{SyllableVector{T}}(undef, 2, 2)

  powers2 = [0, 2]

  rules2[2,2] = id

  rules2[1,2] = [(1, 2), (2, 1)]
  rules2[2,1] = [(1, -2), (2, 1)]

  coll2 = Collector{Int}(findfirst_uncollected_leftmost, powers2, rules2)

  @test normalform(coll2, [(1, 1), (2, 1), (1, 1)]) == [(2, 1)]
  @test normalform(coll2, [(1, 1), (2, 1), (1, -1)]) == [(1, 2), (2, 1)]
  @test normalform(coll2, [(1, 2), (2, 1), (1, -2)]) == [(1, 4), (2, 1)]
end

@testset "SyllableVector{T}" begin
  v = [1, -2, 3, -4]
  sv = SyllableVector{Int}(v)
  @test sv == [(1, 1), (2, -2), (3, 3), (4, -4)]
  @test Vector{Int}(sv, 4) == v

  G = small_group(12, 3)
  id = one(G)
  pcgs = GAP.Globals.Pcgs(G.X)
  x = rand(G)
  y = gens(G)
  v = SyllableVector{Int}(pcgs, x.X)
  @test x == prod([y[i]^e for (i, e) in v]; init = id)

  GG = isomorphic_fp_group(G)[1]
  x = rand(GG)
  y = gens(GG)
  v = SyllableVector{Int}(x)
  @test x == prod([y[i]^e for (i, e) in v])
end

@testset "inv" begin
  @test inv(SyllableVector{Int}([1,2,-1])) == SyllableVector{Int}([(3,1),(2,-2),(1,-1)])
end

@testset "abelianized" begin
  @test abelianized(SyllableVector{Int}()) == []
  @test abelianized([(1, 2), (2, 1), (1, -2)]) == [(2, 1)]
  @test abelianized([(1, 2), (2, -1)]) == [(1, 2), (2, -1)]
  @test abelianized([(2, -1), (1, 2)]) == [(1, 2), (2, -1)]
  @test abelianized([(1, 2), (2, -1), (1, -2), (2, 1)]) == []
end

@testset "collection for random elements in small pc groups" begin
  @test test_collection(Int64, small_group(16, 5))
  @test test_collection(Int32, small_group(32, 8))
  @test test_collection(Int32, small_group(96, 8))
  @test test_collection(Int8, small_group(256, 10000))
end
