module VTLInputs

using FileIO
using GeometryBasics
using LinearAlgebra
using MeshIO
using StaticArrays

export load_stl, compute_sdf, compute_sdf_on_grid, compute_stl_sdf, Triangle, Mesh

"""
    Triangle{T}

A triangle in 3D space defined by three vertices.
"""
struct Triangle{T<:Real}
    v1::SVector{3,T}
    v2::SVector{3,T}
    v3::SVector{3,T}
end

"""
    Mesh{T}

A mesh composed of triangles, loaded from an STL file.
"""
struct Mesh{T<:Real}
    triangles::Vector{Triangle{T}}
end

"""
    load_stl(filepath::AbstractString) -> Mesh{Float64}

Load an STL file and return a Mesh containing all triangles.

# Arguments
- `filepath::AbstractString`: Path to the STL file.

# Returns
- `Mesh{Float64}`: A mesh containing all triangles from the STL file.

# Example
```julia
mesh = load_stl("model.stl")
```
"""
function load_stl(filepath::AbstractString)
    geom = load(filepath)
    
    triangles = Triangle{Float64}[]
    
    # Handle both mesh types from GeometryBasics
    positions = coordinates(geom)
    faces_list = faces(geom)
    
    for face in faces_list
        # Get vertex indices (1-based in Julia)
        idx1, idx2, idx3 = face
        v1 = SVector{3,Float64}(positions[idx1]...)
        v2 = SVector{3,Float64}(positions[idx2]...)
        v3 = SVector{3,Float64}(positions[idx3]...)
        push!(triangles, Triangle(v1, v2, v3))
    end
    
    return Mesh(triangles)
end

"""
    point_to_triangle_distance(point, triangle) -> Float64

Compute the unsigned distance from a point to a triangle in 3D space.
Handles degenerate triangles (zero-area triangles) by computing distance to the nearest vertex.
"""
function point_to_triangle_distance(point::SVector{3,T}, tri::Triangle{T}) where {T<:Real}
    # Edge vectors
    e0 = tri.v2 - tri.v1
    e1 = tri.v3 - tri.v1
    
    # Vector from triangle vertex to point
    v0 = tri.v1 - point
    
    # Compute dot products
    a = dot(e0, e0)
    b = dot(e0, e1)
    c = dot(e1, e1)
    d = dot(e0, v0)
    e = dot(e1, v0)
    
    eps_val = T(1e-12)
    det = a * c - b * b
    
    # Handle degenerate triangles (zero area)
    if abs(det) < eps_val
        # Triangle is degenerate, compute distance to vertices
        d1 = norm(point - tri.v1)
        d2 = norm(point - tri.v2)
        d3 = norm(point - tri.v3)
        return min(d1, d2, d3)
    end
    
    s = b * e - c * d
    t = b * d - a * e
    
    # Safe division helper
    safe_div(num, den) = abs(den) < eps_val ? zero(T) : num / den
    
    if s + t <= det
        if s < 0
            if t < 0
                # Region 4
                if d < 0
                    t = zero(T)
                    s = clamp(safe_div(-d, a), zero(T), one(T))
                else
                    s = zero(T)
                    t = clamp(safe_div(-e, c), zero(T), one(T))
                end
            else
                # Region 3
                s = zero(T)
                t = clamp(safe_div(-e, c), zero(T), one(T))
            end
        elseif t < 0
            # Region 5
            t = zero(T)
            s = clamp(safe_div(-d, a), zero(T), one(T))
        else
            # Region 0 (inside triangle)
            inv_det = one(T) / det
            s *= inv_det
            t *= inv_det
        end
    else
        if s < 0
            # Region 2
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0
                numer = tmp1 - tmp0
                denom = a - 2 * b + c
                s = clamp(safe_div(numer, denom), zero(T), one(T))
                t = one(T) - s
            else
                s = zero(T)
                t = clamp(safe_div(-e, c), zero(T), one(T))
            end
        elseif t < 0
            # Region 6
            tmp0 = b + e
            tmp1 = a + d
            if tmp1 > tmp0
                numer = tmp1 - tmp0
                denom = a - 2 * b + c
                t = clamp(safe_div(numer, denom), zero(T), one(T))
                s = one(T) - t
            else
                t = zero(T)
                s = clamp(safe_div(-d, a), zero(T), one(T))
            end
        else
            # Region 1
            numer = (c + e) - (b + d)
            denom = a - 2 * b + c
            s = clamp(safe_div(numer, denom), zero(T), one(T))
            t = one(T) - s
        end
    end
    
    # Closest point on triangle
    closest = tri.v1 + s * e0 + t * e1
    
    return norm(point - closest)
end

"""
    triangle_normal(triangle) -> SVector{3}

Compute the outward normal of a triangle.
"""
function triangle_normal(tri::Triangle{T}) where {T<:Real}
    e0 = tri.v2 - tri.v1
    e1 = tri.v3 - tri.v1
    n = cross(e0, e1)
    len = norm(n)
    return len > eps(T) ? n / len : SVector{3,T}(0, 0, 1)
end

"""
    ray_triangle_intersection(origin, direction, triangle) -> (Bool, Float64)

Check if a ray intersects a triangle using the Möller–Trumbore algorithm.
Returns (intersects, t) where t is the distance along the ray.
"""
function ray_triangle_intersection(origin::SVector{3,T}, direction::SVector{3,T}, tri::Triangle{T}) where {T<:Real}
    eps_val = T(1e-10)
    
    e1 = tri.v2 - tri.v1
    e2 = tri.v3 - tri.v1
    h = cross(direction, e2)
    a = dot(e1, h)
    
    if abs(a) < eps_val
        return (false, T(0))
    end
    
    f = one(T) / a
    s = origin - tri.v1
    u = f * dot(s, h)
    
    if u < 0 || u > 1
        return (false, T(0))
    end
    
    q = cross(s, e1)
    v = f * dot(direction, q)
    
    if v < 0 || u + v > 1
        return (false, T(0))
    end
    
    t = f * dot(e2, q)
    
    if t > eps_val
        return (true, t)
    end
    
    return (false, T(0))
end

"""
    compute_sign_raycast(point, mesh) -> Int

Compute the sign of the distance using ray casting.
Returns -1 if point is inside the mesh, +1 if outside.
Uses multiple rays for robustness.
"""
function compute_sign_raycast(point::SVector{3,T}, mesh::Mesh{T}) where {T<:Real}
    # Cast rays in multiple directions for robustness
    directions = [
        SVector{3,T}(1, 0, 0),
        SVector{3,T}(0, 1, 0),
        SVector{3,T}(0, 0, 1),
        SVector{3,T}(1, 1, 1) / sqrt(T(3)),
        SVector{3,T}(-1, 1, 1) / sqrt(T(3))
    ]
    
    votes_inside = 0
    votes_outside = 0
    eps_val = T(1e-8)
    
    for dir in directions
        # Collect unique intersection distances (merge co-planar triangles)
        intersection_ts = T[]
        for tri in mesh.triangles
            intersects, t = ray_triangle_intersection(point, dir, tri)
            if intersects && t > eps_val
                # Check if this is a new unique intersection
                is_unique = true
                for existing_t in intersection_ts
                    if abs(t - existing_t) < eps_val
                        is_unique = false
                        break
                    end
                end
                if is_unique
                    push!(intersection_ts, t)
                end
            end
        end
        
        # Odd number of unique intersections means inside
        if isodd(length(intersection_ts))
            votes_inside += 1
        else
            votes_outside += 1
        end
    end
    
    return votes_inside > votes_outside ? -1 : 1
end

"""
    compute_sign(point, triangle) -> Int

Compute the sign of the distance based on the triangle normal.
Returns -1 if point is inside (behind the triangle), +1 if outside.
"""
function compute_sign(point::SVector{3,T}, tri::Triangle{T}) where {T<:Real}
    n = triangle_normal(tri)
    # Vector from triangle centroid to point
    centroid = (tri.v1 + tri.v2 + tri.v3) / 3
    to_point = point - centroid
    return dot(n, to_point) >= 0 ? 1 : -1
end

"""
    compute_sdf(point, mesh) -> Float64

Compute the signed distance from a point to a mesh.
Negative values indicate the point is inside the mesh,
positive values indicate the point is outside.

# Arguments
- `point`: A 3D point as a tuple, array, or SVector.
- `mesh::Mesh`: The mesh to compute distance to.

# Returns
- `Float64`: The signed distance from the point to the mesh.

# Example
```julia
mesh = load_stl("model.stl")
sdf = compute_sdf((0.0, 0.0, 0.0), mesh)
```
"""
function compute_sdf(point, mesh::Mesh{T}) where {T<:Real}
    p = SVector{3,T}(point...)
    
    # Handle empty mesh
    if isempty(mesh.triangles)
        return T(Inf)
    end
    
    min_dist = T(Inf)
    closest_tri_idx = 1
    
    # Find closest triangle
    for (i, tri) in enumerate(mesh.triangles)
        dist = point_to_triangle_distance(p, tri)
        if dist < min_dist
            min_dist = dist
            closest_tri_idx = i
        end
    end
    
    # Compute sign
    # Use ray casting for meshes with multiple triangles (closed surfaces)
    # Use normal-based sign for single triangles
    if length(mesh.triangles) > 1
        sign = compute_sign_raycast(p, mesh)
    else
        sign = compute_sign(p, mesh.triangles[closest_tri_idx])
    end
    
    return sign * min_dist
end

"""
    compute_sdf_on_grid(grid_points, mesh) -> Array{Float64}

Compute the signed distance function on a grid of points.
Works with NTuple vectors representing a background mesh.

# Arguments
- `grid_points`: An iterable of 3D points (can be tuples, arrays, or SVectors).
- `mesh::Mesh`: The mesh to compute distances to.

# Returns
- `Vector{Float64}`: A vector of signed distances for each grid point.

# Example
```julia
mesh = load_stl("model.stl")
# Create a grid of points
xs = range(-1, 1, length=10)
ys = range(-1, 1, length=10)
zs = range(-1, 1, length=10)
grid = [(x, y, z) for x in xs for y in ys for z in zs]
sdf_values = compute_sdf_on_grid(grid, mesh)
```
"""
function compute_sdf_on_grid(grid_points, mesh::Mesh{T}) where {T<:Real}
    return [compute_sdf(p, mesh) for p in grid_points]
end

"""
    compute_sdf_on_grid(xs, ys, zs, mesh) -> Array{Float64, 3}

Compute the signed distance function on a 3D structured grid.

# Arguments
- `xs`: X coordinates of the grid.
- `ys`: Y coordinates of the grid.
- `zs`: Z coordinates of the grid.
- `mesh::Mesh`: The mesh to compute distances to.

# Returns
- `Array{Float64, 3}`: A 3D array of signed distances.

# Example
```julia
mesh = load_stl("model.stl")
xs = range(-1, 1, length=10)
ys = range(-1, 1, length=10)
zs = range(-1, 1, length=10)
sdf_grid = compute_sdf_on_grid(xs, ys, zs, mesh)
```
"""
function compute_sdf_on_grid(xs, ys, zs, mesh::Mesh{T}) where {T<:Real}
    nx, ny, nz = length(xs), length(ys), length(zs)
    sdf = Array{T}(undef, nx, ny, nz)
    
    for (i, x) in enumerate(xs)
        for (j, y) in enumerate(ys)
            for (k, z) in enumerate(zs)
                sdf[i, j, k] = compute_sdf((x, y, z), mesh)
            end
        end
    end
    
    return sdf
end

"""
    compute_stl_sdf(filepath::AbstractString; dims::Int=3) -> Function

Load an STL file and return a signed distance function.

# Arguments
- `filepath::AbstractString`: Path to the STL file.
- `dims::Int=3`: Dimensionality of the SDF function. Use `3` for 3D `(x, y, z) -> sdf` 
  or `2` for 2D `(x, y) -> sdf` (projects to z=0 plane).

# Returns
- For `dims=3`: A function `(x, y, z) -> Float64` that computes the signed distance.
- For `dims=2`: A function `(x, y) -> Float64` that computes the signed distance at z=0.

# Example
```julia
# 3D SDF function
sdf_3d = compute_stl_sdf("model.stl")
distance = sdf_3d(0.5, 0.5, 0.5)

# 2D SDF function (evaluates at z=0)
sdf_2d = compute_stl_sdf("model.stl", dims=2)
distance = sdf_2d(0.5, 0.5)
```
"""
function compute_stl_sdf(filepath::AbstractString; dims::Int=3)
    mesh = load_stl(filepath)
    return compute_stl_sdf(mesh; dims=dims)
end

"""
    compute_stl_sdf(mesh::Mesh; dims::Int=3) -> Function

Create a signed distance function from a mesh.

# Arguments
- `mesh::Mesh`: The mesh to compute distances to.
- `dims::Int=3`: Dimensionality of the SDF function. Use `3` for 3D `(x, y, z) -> sdf` 
  or `2` for 2D `(x, y) -> sdf` (projects to z=0 plane).

# Returns
- For `dims=3`: A function `(x, y, z) -> Float64` that computes the signed distance.
- For `dims=2`: A function `(x, y) -> Float64` that computes the signed distance at z=0.

# Example
```julia
mesh = load_stl("model.stl")

# 3D SDF function
sdf_3d = compute_stl_sdf(mesh)
distance = sdf_3d(0.5, 0.5, 0.5)

# 2D SDF function (evaluates at z=0)
sdf_2d = compute_stl_sdf(mesh, dims=2)
distance = sdf_2d(0.5, 0.5)
```
"""
function compute_stl_sdf(mesh::Mesh{T}; dims::Int=3) where {T<:Real}
    if dims == 3
        return (x, y, z) -> compute_sdf((x, y, z), mesh)
    elseif dims == 2
        return (x, y) -> compute_sdf((x, y, T(0)), mesh)
    else
        throw(ArgumentError("dims must be 2 or 3, got $dims"))
    end
end

end # module
