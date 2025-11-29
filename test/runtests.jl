using Test
using VTLInputs
using StaticArrays
using LinearAlgebra

# Get the path to test data directory
const TESTDATA_DIR = joinpath(@__DIR__, "testdata")

@testset "VTLInputs.jl" begin
    
    @testset "Triangle creation" begin
        v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v3 = SVector{3,Float64}(0.0, 1.0, 0.0)
        tri = Triangle(v1, v2, v3)
        
        @test tri.v1 == v1
        @test tri.v2 == v2
        @test tri.v3 == v3
    end
    
    @testset "Mesh creation" begin
        v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v3 = SVector{3,Float64}(0.0, 1.0, 0.0)
        tri = Triangle(v1, v2, v3)
        mesh = Mesh([tri])
        
        @test length(mesh.triangles) == 1
        @test mesh.triangles[1] == tri
    end
    
    @testset "SDF computation for simple triangle" begin
        # Create a simple triangle in the XY plane
        v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v3 = SVector{3,Float64}(0.5, 1.0, 0.0)
        tri = Triangle(v1, v2, v3)
        mesh = Mesh([tri])
        
        # Point above the triangle (positive z)
        point_above = (0.5, 0.5, 1.0)
        sdf_above = compute_sdf(point_above, mesh)
        @test sdf_above > 0  # Should be positive (outside)
        
        # Point below the triangle (negative z)
        point_below = (0.5, 0.5, -1.0)
        sdf_below = compute_sdf(point_below, mesh)
        @test sdf_below < 0  # Should be negative (inside)
        
        # Point on the triangle
        point_on = (0.5, 0.5, 0.0)
        sdf_on = compute_sdf(point_on, mesh)
        @test abs(sdf_on) < 0.1  # Should be near zero
    end
    
    @testset "SDF on grid" begin
        # Create a simple triangle mesh
        v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v3 = SVector{3,Float64}(0.5, 1.0, 0.0)
        tri = Triangle(v1, v2, v3)
        mesh = Mesh([tri])
        
        # Test with tuple-based grid points
        grid_points = [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.5, 0.5, 0.0)]
        sdf_values = compute_sdf_on_grid(grid_points, mesh)
        
        @test length(sdf_values) == 3
        @test sdf_values[1] > 0  # Above triangle
        @test sdf_values[2] < 0  # Below triangle
        @test abs(sdf_values[3]) < 0.1  # On triangle
    end
    
    @testset "SDF on 3D structured grid" begin
        # Create a simple triangle mesh
        v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(1.0, 0.0, 0.0)
        v3 = SVector{3,Float64}(0.5, 1.0, 0.0)
        tri = Triangle(v1, v2, v3)
        mesh = Mesh([tri])
        
        # Create a small 3D grid
        xs = range(-1, 2, length=3)
        ys = range(-1, 2, length=3)
        zs = range(-1, 1, length=3)
        
        sdf_grid = compute_sdf_on_grid(xs, ys, zs, mesh)
        
        @test size(sdf_grid) == (3, 3, 3)
    end
    
    @testset "Distance computation accuracy" begin
        # Create a triangle at z=0
        v1 = SVector{3,Float64}(0.0, 0.0, 0.0)
        v2 = SVector{3,Float64}(2.0, 0.0, 0.0)
        v3 = SVector{3,Float64}(1.0, 2.0, 0.0)
        tri = Triangle(v1, v2, v3)
        mesh = Mesh([tri])
        
        # Point directly above the centroid
        centroid_x = (0.0 + 2.0 + 1.0) / 3
        centroid_y = (0.0 + 0.0 + 2.0) / 3
        height = 5.0
        point = (centroid_x, centroid_y, height)
        
        sdf = compute_sdf(point, mesh)
        @test abs(abs(sdf) - height) < 0.1  # Distance should be approximately the height
    end
    
    @testset "STL file loading" begin
        cube_path = joinpath(TESTDATA_DIR, "cube.stl")
        mesh = load_stl(cube_path)
        
        # A cube has 12 triangles (2 per face, 6 faces)
        @test length(mesh.triangles) == 12
        
        # Test SDF for the cube
        # Point at center of cube should be inside (negative)
        center_point = (0.5, 0.5, 0.5)
        sdf_center = compute_sdf(center_point, mesh)
        @test sdf_center < 0
        
        # Point far outside cube should be outside (positive)
        outside_point = (5.0, 5.0, 5.0)
        sdf_outside = compute_sdf(outside_point, mesh)
        @test sdf_outside > 0
    end
    
    @testset "Edge cases" begin
        # Empty mesh
        empty_mesh = Mesh(Triangle{Float64}[])
        @test compute_sdf((0.0, 0.0, 0.0), empty_mesh) == Inf
        
        # Degenerate triangle (zero area - all vertices same)
        v = SVector{3,Float64}(1.0, 1.0, 1.0)
        degenerate_tri = Triangle(v, v, v)
        degenerate_mesh = Mesh([degenerate_tri])
        sdf = compute_sdf((0.0, 0.0, 0.0), degenerate_mesh)
        @test isfinite(sdf)  # Should not crash
        @test abs(abs(sdf) - sqrt(3.0)) < 0.01  # Distance from origin to (1,1,1) should be sqrt(3)
    end
    
end
