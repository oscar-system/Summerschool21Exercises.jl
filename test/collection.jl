@testset "entered by hand: D8" begin
  id = SyllableVector{Int}()
  rules = Base.Matrix{SyllableVector{Int}}(undef, 3, 3)

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
  
  coll = Collector{Int}(findfirst_uncollected_leftmost, powers, rules)
  
  @test normalform(coll, [(1, 2)]) == [(2, 1)]
  @test normalform(coll, [(1, 4)]) == id
  @test normalform(coll, [(1, 1), (2, 1), (1, 1)]) == id
  @test normalform(coll, [(3, 1), (1, 1)]) == [(1, 1), (2, 1), (3, 1)]
end

@testset "SyllableVector{T}" begin
  v = [1, -2, 3, -4]
  sv = SyllableVector{Int}(v)
  @test sv == [(1, 1), (2, -2), (3, 3), (4, -4)]
  @test Vector{Int}(sv, 4) == v

  G = small_group(12, 3)
  pcgs = GAP.Globals.Pcgs(G.X)
  x = rand(G)
  y = gens(G)
  v = SyllableVector{Int}(pcgs, x.X)
  @test x == prod([y[i]^e for (i, e) in v])

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
