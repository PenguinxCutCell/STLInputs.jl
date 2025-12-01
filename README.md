# VTLInputs.jl

[![Build Status](https://github.com/PenguinxCutCell/VTLInputs.jl/workflows/CI/badge.svg)](https://github.com/PenguinxCutCell/VTLInputs.jl/actions)

A Julia package for reading STL files and computing signed distance functions (SDF) for 3D objects.

## Installation

```julia
using Pkg
Pkg.add(url="https://github.com/PenguinxCutCell/VTLInputs.jl")
```

## Features

- Load STL files (ASCII and binary formats)
- Compute signed distance function (SDF) for points relative to mesh surfaces
- Compute SDF on grids (both arbitrary point collections and structured 3D grids)
- **Create SDF function objects** for easy evaluation at any point
- **Support for both 2D and 3D** SDF computations
- Negative SDF values indicate points inside the mesh, positive values indicate outside

## Quick Start

```julia
using VTLInputs

# Load an STL file and get an SDF function object
sdf = compute_stl_sdf("model.stl")

# Evaluate SDF at any point using the function
distance = sdf(0.5, 0.5, 0.5)  # 3D evaluation

# For 2D applications (evaluates at z=0)
sdf_2d = compute_stl_sdf("model.stl", dims=2)
distance_2d = sdf_2d(0.5, 0.5)  # 2D evaluation
```

### Alternative: Direct Computation

```julia
using VTLInputs

# Load an STL file
mesh = load_stl("model.stl")

# Compute SDF at a single point
point = (0.5, 0.5, 0.5)
sdf_value = compute_sdf(point, mesh)

# Compute SDF on a collection of points (NTuple vectors from a background mesh)
grid_points = [(x, y, z) for x in 0:0.1:1 for y in 0:0.1:1 for z in 0:0.1:1]
sdf_values = compute_sdf_on_grid(grid_points, mesh)

# Compute SDF on a structured 3D grid
xs = range(-1, 1, length=50)
ys = range(-1, 1, length=50)
zs = range(-1, 1, length=50)
sdf_grid = compute_sdf_on_grid(xs, ys, zs, mesh)
```

## API Reference

### Types

- `Triangle{T}`: A triangle in 3D space defined by three vertices (v1, v2, v3)
- `Mesh{T}`: A mesh composed of triangles, loaded from an STL file

### Functions

- `load_stl(filepath)`: Load an STL file and return a Mesh
- `compute_stl_sdf(filepath; dims=3)`: Load an STL file and return an SDF function object
  - `dims=3`: Returns `(x, y, z) -> sdf_value`
  - `dims=2`: Returns `(x, y) -> sdf_value` (evaluates at z=0)
- `compute_stl_sdf(mesh; dims=3)`: Create an SDF function object from a mesh
- `compute_sdf(point, mesh)`: Compute the signed distance from a point to a mesh
- `compute_sdf_on_grid(grid_points, mesh)`: Compute SDF for a collection of 3D points
- `compute_sdf_on_grid(xs, ys, zs, mesh)`: Compute SDF on a structured 3D grid

## Usage with Background Meshes

This package is designed to work with background meshes represented as NTuple vectors:

```julia
using VTLInputs

# Get an SDF function object
sdf = compute_stl_sdf("geometry.stl")

# Evaluate at any point
inside = sdf(1.0, 2.0, 3.0) < 0  # true if point is inside

# Or use the grid-based approach
mesh = load_stl("geometry.stl")

# Define your background mesh as NTuple vectors
nx, ny, nz = 100, 100, 100
xs = range(-5, 5, length=nx)
ys = range(-5, 5, length=ny)
zs = range(-5, 5, length=nz)

# Compute SDF on the background mesh
sdf_grid = compute_sdf_on_grid(xs, ys, zs, mesh)

# Now sdf_grid[i,j,k] contains the signed distance at position (xs[i], ys[j], zs[k])
# Negative values: inside the object
# Positive values: outside the object
# Zero (or near zero): on the surface
```

## 2D SDF Support

For 2D applications or when working with planar STL files:

```julia
using VTLInputs

# Create a 2D SDF function (evaluates at z=0)
sdf_2d = compute_stl_sdf("model.stl", dims=2)

# Evaluate at 2D points
distance = sdf_2d(x, y)
```

## License

MIT