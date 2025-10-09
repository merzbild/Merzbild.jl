@muladd begin

"""
    SurfProps

Structure to store computed surface properties.

# Fields
* `n_elements`: number of physical cells
* `n_species`: number of species
* `areas`: the vector of surface element areas
* `inv_areas`: the vector of the inverse surface element areas
* `normals`: the array of surface element normals with shape  `(3, n_elements)`
* `np`: number of particles impacting the surface elements, array of shape `(n_elements, n_species)`
* `flux_incident`: incident mass flux per surface element, array of shape `(n_elements, n_species)`
* `flux_reflected`: reflected mass flux per surface element, array of shape `(n_elements, n_species)`
* `force`: force acting per surface element, array of shape `(3, n_elements, n_species)`
* `normal_pressure`: normal pressure per surface element, array of shape `(n_elements, n_species)`
* `shear_pressure`: shear pressure per surface element, array of shape `(3, n_elements, n_species)`
* `kinetic_energy_flux`: kinetic energy flux per surface element, array of shape `(n_elements, n_species)`
"""
mutable struct SurfProps
    n_elements::Int64
    n_species::Int64
    areas::Array{Float64}  # surface element areas: elements
    inv_areas::Array{Float64}  # surface element areas: elements
    normals::Array{Float64,2}  # surface normals: component x elements
    np::Array{Float64,2}  # number of particles: elements x species
    flux_incident::Array{Float64,2}  # incident mass flux: elements x species
    flux_reflected::Array{Float64,2}  # reflected mass flux: elements x species
    force::Array{Float64,3}  # force: force component x elements x species
    normal_pressure::Array{Float64,2}  # pressure: elements x species
    shear_pressure::Array{Float64,3}  # shear pressure: shear pressure component x elements x species
    kinetic_energy_flux::Array{Float64,2}  # kinetic energy flux: elements x species

    @doc """
        SurfProps(n_elements, n_species, area, normals)

    # Positional arguments
    * `n_elements`
    * `n_species`
    * `areas`
    * `normals`
    """
    function SurfProps(n_elements, n_species, areas, normals)
        return new(n_elements, n_species, areas, [1.0/area for area in areas], normals,
                   zeros(n_elements, n_species), zeros(n_elements, n_species), zeros(n_elements, n_species),
                   zeros(3, n_elements, n_species), zeros(n_elements, n_species), zeros(3, n_elements, n_species), zeros(n_elements, n_species))
    end
end

"""
    SurfProps(pia, grid::Grid1DUniform)

Create a `SurfProps` struct for a 1-D grid, with element 1 corresponding to the left wall
and element 2 corresponding to the right wall. The areas of the wall are assumed to
be equal to 1, the normals are parallel to the x axis.

# Positional arguments
* `pia`: the `ParticleIndexerArray` instance
* `grid`: the `Grid1DUniform` grid
"""
SurfProps(pia, grid::Grid1DUniform) = SurfProps(2, size(pia.indexer)[2], [1.0, 1.0],
                                                      [1.0 0.0 0.0; -1.0 0.0 0.0]')

"""
    update_surface_incident!(particle, species, surf_props, surface_element_id)

Update surface properties for surface element `surface_element_id` for an incident particle.

# Positional arguments
* `particle`: the particle hitting the surface before its velocity is updated
* `species`: the species of the particle
* `surf_props`: the `SurfProps` instance
* `surface_element_id`: id of the surface element which the particle has impacted
"""
function update_surface_incident!(particle, species, surf_props, surface_element_id)
    @inbounds px = particle.w * particle.v[1]
    @inbounds py = particle.w * particle.v[2]
    @inbounds pz = particle.w * particle.v[3]

    @inbounds p_dot_n = px * surf_props.normals[1, surface_element_id] + py * surf_props.normals[2, surface_element_id] + pz * surf_props.normals[3, surface_element_id]

    @inbounds surf_props.np[surface_element_id, species] += 1
    @inbounds surf_props.flux_incident[surface_element_id, species] += particle.w

    @inbounds surf_props.force[1, surface_element_id, species] += px
    @inbounds surf_props.force[2, surface_element_id, species] += py
    @inbounds surf_props.force[3, surface_element_id, species] += pz

    @inbounds surf_props.normal_pressure[surface_element_id, species] -= p_dot_n

    @inbounds surf_props.shear_pressure[1, surface_element_id, species] += (px - p_dot_n * surf_props.normals[1, surface_element_id])
    @inbounds surf_props.shear_pressure[2, surface_element_id, species] += (py - p_dot_n * surf_props.normals[2, surface_element_id])
    @inbounds surf_props.shear_pressure[3, surface_element_id, species] += (pz - p_dot_n * surf_props.normals[3, surface_element_id])

    @inbounds surf_props.kinetic_energy_flux[surface_element_id, species] += 0.5 * (px * particle.v[1] + py * particle.v[2] + pz * particle.v[3])
end

"""
    update_surface_incident!(particle, species, surf_props, surface_element_id)

Update surface properties for surface element `surface_element_id` for a reflected particle.

# Positional arguments
* `particle`: the particle hitting the surface after its velocity is updated
* `species`: the species of the particle
* `surf_props`: the `SurfProps` instance
* `surface_element_id`: id of the surface element which the particle has impacted
"""
function update_surface_reflected!(particle, species, surf_props, surface_element_id)

    @inbounds px = particle.w * particle.v[1]
    @inbounds py = particle.w * particle.v[2]
    @inbounds pz = particle.w * particle.v[3]

    @inbounds p_dot_n = px * surf_props.normals[1, surface_element_id] + py * surf_props.normals[2, surface_element_id] + pz * surf_props.normals[3, surface_element_id]

    @inbounds surf_props.flux_reflected[surface_element_id, species] -= particle.w

    @inbounds surf_props.force[1, surface_element_id, species] -= px
    @inbounds surf_props.force[2, surface_element_id, species] -= py
    @inbounds surf_props.force[3, surface_element_id, species] -= pz

    @inbounds surf_props.normal_pressure[surface_element_id, species] += p_dot_n
    @inbounds surf_props.shear_pressure[1, surface_element_id, species] -= (px - p_dot_n * surf_props.normals[1, surface_element_id])
    @inbounds surf_props.shear_pressure[2, surface_element_id, species] -= (py - p_dot_n * surf_props.normals[2, surface_element_id])
    @inbounds surf_props.shear_pressure[3, surface_element_id, species] -= (pz - p_dot_n * surf_props.normals[3, surface_element_id])

    @inbounds surf_props.kinetic_energy_flux[surface_element_id, species] -= 0.5 * (px * particle.v[1] + py * particle.v[2] + pz * particle.v[3])
end

"""
    surface_props_scale!(species, surf_props, species_data, Δt)

Scale computed surface properties using the molecular mass of species, the inverse of the timestep, and the inverse surface area.

# Positional arguments
* `species`: the species of the particle
* `surf_props`: the `SurfProps` instance
* `species_data`: the vector of `SpeciesData` of the chemical species in the flow
* `Δt`: the timestep over which the surface properties were computed
"""
function surface_props_scale!(species, surf_props, species_data, Δt)
    @inbounds factor_base = species_data[species].mass / Δt

    @inbounds @simd for surface_element_id in 1:surf_props.n_elements
        factor = factor_base * surf_props.inv_areas[surface_element_id]

        surf_props.flux_incident[surface_element_id, species] *= factor
        surf_props.flux_reflected[surface_element_id, species] *= factor

        for i in 1:3
            surf_props.force[i,surface_element_id,species] *= factor
            surf_props.shear_pressure[i,surface_element_id,species] *= factor
        end
        surf_props.normal_pressure[surface_element_id,species] *= factor
        surf_props.kinetic_energy_flux[surface_element_id,species] *= factor
    end
end

"""
    clear_props!(surf_props::SurfProps)

Clear all data from a `SurfProps` instance, either at the start of a new convection step, or when physical properties are averaged over timesteps
and averaging over a new set of timesteps needs to be started.

# Positional arguments
* `surf_props`: the `SurfProps` instance to be cleared
"""
function clear_props!(surf_props::SurfProps)
    surf_props.np[:,:] .= 0
    surf_props.flux_incident[:,:] .= 0.0
    surf_props.flux_reflected[:,:] .= 0.0
    surf_props.force[:,:,:] .= 0.0
    surf_props.normal_pressure[:,:] .= 0.0
    surf_props.shear_pressure[:,:,:] .= 0.0
    surf_props.kinetic_energy_flux[:,:] .= 0.0
end

"""
    avg_props!(surf_props_avg::SurfProps, surf_props::SurfProps, n_avg_timesteps)

Used to time-average computed surface properties.
For each instantaneous value of a property computed and stored in `surf_props`,
it is divided by `n_avg_timesteps` and added to `surf_props_avg`.

# Positional arguments
* `surf_props_avg`: the `SurfProps` instance used to store the time-averaged properties
* `surf_props`: the `SurfProps` instance holding the current values of the properties
    to be used for the averaging at the current timestep
* `n_avg_timesteps`: the number of timesteps over which the averaging is performed
"""
function avg_props!(surf_props_avg::SurfProps, surf_props::SurfProps, n_avg_timesteps)
    inv_nt_avg = 1.0 / n_avg_timesteps

    for species in 1:surf_props.n_species
        @inbounds @simd for element in 1:surf_props.n_elements
            surf_props_avg.np[element,species] = surf_props_avg.np[element,species] + surf_props.np[element,species] * inv_nt_avg
            surf_props_avg.flux_incident[element,species] = surf_props_avg.flux_incident[element,species] + surf_props.flux_incident[element,species] * inv_nt_avg
            surf_props_avg.flux_reflected[element,species] = surf_props_avg.flux_reflected[element,species] + surf_props.flux_reflected[element,species] * inv_nt_avg
            surf_props_avg.force[1,element,species] = surf_props_avg.force[1,element,species] + surf_props.force[1,element,species] * inv_nt_avg
            surf_props_avg.force[2,element,species] = surf_props_avg.force[2,element,species] + surf_props.force[2,element,species] * inv_nt_avg
            surf_props_avg.force[3,element,species] = surf_props_avg.force[3,element,species] + surf_props.force[3,element,species] * inv_nt_avg
            surf_props_avg.normal_pressure[element,species] = surf_props_avg.normal_pressure[element,species] + surf_props.normal_pressure[element,species] * inv_nt_avg
            surf_props_avg.shear_pressure[1,element,species] = surf_props_avg.shear_pressure[1,element,species] + surf_props.shear_pressure[1,element,species] * inv_nt_avg
            surf_props_avg.shear_pressure[2,element,species] = surf_props_avg.shear_pressure[2,element,species] + surf_props.shear_pressure[2,element,species] * inv_nt_avg
            surf_props_avg.shear_pressure[3,element,species] = surf_props_avg.shear_pressure[3,element,species] + surf_props.shear_pressure[3,element,species] * inv_nt_avg
            surf_props_avg.kinetic_energy_flux[element,species] = surf_props_avg.kinetic_energy_flux[element,species] + surf_props.kinetic_energy_flux[element,species] * inv_nt_avg
        end
    end
end

"""
    reduce_surf_props!(surf_props_target, surf_props_chunks)

Sum up the values of the computed surface properties for all `SurfProps` instances in 
a `surf_props_chunks` list and store the sums in `surf_props_target`.

# Positional arguments
* `surf_props_target`: the `SurfProps` instance which will hold the reduced values
* `surf_props_chunks`: the list `SurfProps` instances to use for the reduction operation
"""
function reduce_surf_props!(surf_props_target, surf_props_chunks)
    clear_props!(surf_props_target)

    for surf_props in surf_props_chunks
        for species in 1:surf_props.n_species
            @inbounds @simd for element in 1:surf_props.n_elements
                surf_props_target.np[element,species] += surf_props.np[element,species]
                surf_props_target.flux_incident[element,species] += surf_props.flux_incident[element,species]
                surf_props_target.flux_reflected[element,species] += surf_props.flux_reflected[element,species]
                surf_props_target.force[1,element,species] += surf_props.force[1,element,species]
                surf_props_target.force[2,element,species] += surf_props.force[2,element,species]
                surf_props_target.force[3,element,species] += surf_props.force[3,element,species]
                surf_props_target.normal_pressure[element,species] += surf_props.normal_pressure[element,species]
                surf_props_target.shear_pressure[1,element,species] += surf_props.shear_pressure[1,element,species]
                surf_props_target.shear_pressure[2,element,species] += surf_props.shear_pressure[2,element,species]
                surf_props_target.shear_pressure[3,element,species] += surf_props.shear_pressure[3,element,species]
                surf_props_target.kinetic_energy_flux[element,species] += surf_props.kinetic_energy_flux[element,species]
            end
        end
    end
end

end