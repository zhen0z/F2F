export Patchy

# Define Patchy for 3D arrays (spatial)
function Patchy(sg::AbstractArray{T, 3}, phi::AbstractArray{T, 3};
    bulk_grain::T = T(39.3f9), bulk_frame::T = T(2.56f9), shear_frame::T = T(8.5f9), bulk_fl1::T = T(2.31f9), 
    bulk_fl2::T = T(0.08f9), rhoG::T = T(700.0f0), rhoL::T = T(1190.0f0), rhoR::T = T(2100.0f0)) where T

    bulk_sat1 = bulk_frame .+ ((1 .- bulk_frame ./ bulk_grain) .^ 2) ./ ((phi ./ bulk_fl1) .+ ((1 .- phi) ./ bulk_grain) .- bulk_frame ./ bulk_grain .^ 2)
    rhoR1 = rhoR .* (1 .- phi) .+ phi .* rhoL

    patch_temp = bulk_sat1 ./ (bulk_grain .- bulk_sat1) .- 
                 bulk_fl1 ./ phi ./ (bulk_grain .- bulk_fl1) .+ 
                 bulk_fl2 ./ phi ./ (bulk_grain .- bulk_fl2)

    bulk_sat2 = bulk_grain ./ (1f0 ./ patch_temp .+ 1f0)

    bulk_new = 1f0 ./ ((1f0 .- sg) ./ (bulk_sat1 .+ 4f0 ./ 3f0 .* shear_frame) .+ 
                       sg ./ (bulk_sat2 .+ 4f0 ./ 3f0 .* shear_frame)) .- 4f0/3f0 .* shear_frame

    rho_new = rhoR .+ phi .* sg .* (rhoG .- rhoL)

    Vp_new = sqrt.((bulk_new .+ 4f0 ./ 3f0 .* shear_frame) ./ rho_new)
    return Vp_new, rho_new
end

# Define Patchy for 4D arrays (time steps)
function Patchy(sg::AbstractArray{T, 4}, phi::AbstractArray{T, 3};
    bulk_grain::T = T(39.3f9), bulk_frame::T = T(2.56f9), shear_frame::T = T(8.5f9), bulk_fl1::T = T(2.31f9), 
    bulk_fl2::T = T(0.08f9), rhoG::T = T(700.0f0), rhoL::T = T(1190.0f0), rhoR::T = T(2100.0f0)) where T

    # Loop over the time dimension (1st dimension of sg)
    stack = [Patchy(sg[i, :, :, :], phi, bulk_grain=bulk_grain; bulk_frame=bulk_frame, shear_frame=shear_frame, bulk_fl1=bulk_fl1, bulk_fl2=bulk_fl2, rhoG=rhoG, rhoL=rhoL, rhoR=rhoR) for i = 1:size(sg, 1)]

    # Collect results for each time step
    Vp = [stack[i][1] for i = 1:size(sg, 1)]
    ρ = [stack[i][2] for i = 1:size(sg, 1)]

    return Vp, ρ
end
