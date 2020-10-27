using AbstractPlotting
using MakieLayout
import AbstractPlotting: textslider
using GLMakie

z_slider, z_obs = textslider(LinRange(1.0, 10.0, 101), "z max", start = 9.0)
v_slider, v_obs = textslider(LinRange(0.0, 10.0, 101), "v", start = 1.5)
w_slider, w_obs = textslider(LinRange(0.0, 10.0, 101), "w", start = 2.0)
β_slider, β_obs = textslider(LinRange(0.0, 10.0, 101), "β", start = 8.0)

vx = -20:0.01:20
vy = -20:0.01:20

ℜD = lift(z_obs, v_obs, w_obs, β_obs) do z_max, v, w, β
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    f(x, y) = R * cos(v * x) * (exp(-v * y) + 2 * P * cosh(v * y)) - (x^2 - y^2) / β - y - R * (2 * P + 1)
    [ifelse(abs(f(x, y)) < z_max, f(x, y), NaN) for x = vx, y = vy]
end

ℑD = lift(z_obs, v_obs, w_obs, β_obs) do z_max, v, w, β
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    f(x, y) = R * sin(v * x) * (exp(-v * y) - 2 * P * sinh(v * y)) - x * (2 * y / β - 1)
    [ifelse(abs(f(x, y)) < z_max, f(x, y), NaN) for x = vx, y = vy]
end

vx_mesh = -20:0.25:20
vy_mesh = -20:0.25:20

ℜD_mesh = lift(z_obs, v_obs, w_obs, β_obs) do z_max, v, w, β
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    f(x, y) = R * cos(v * x) * (exp(-v * y) + 2 * P * cosh(v * y)) - (x^2 - y^2) / β - y - R * (2 * P + 1)
    [ifelse(abs(f(x, y)) < z_max, f(x, y), NaN) for x = vx_mesh, y = vy_mesh]
end

ℑD_mesh = lift(z_obs, v_obs, w_obs, β_obs) do z_max, v, w, β
    R = (v^2 - w^2) / (w^2 * v)
    P = 1 / (exp(v * β) - 1)
    f(x, y) = R * sin(v * x) * (exp(-v * y) - 2 * P * sinh(v * y)) - x * (2 * y / β - 1)
    [ifelse(abs(f(x, y)) < z_max, f(x, y), NaN) for x = vx_mesh, y = vy_mesh]
end

scene_real = Scene(resolution = (1000, 1000))
surface!(scene_real, vx, vy, ℜD)
xmr, ymr, zmr = minimum(scene_limits(scene_real))
contour!(scene_real, vx, vy, ℜD, levels = 30, linewidth = 1, transformation = (:xy, zmr))
wireframe!(scene_real, vx_mesh, vy_mesh, ℜD_mesh, overdraw = true, transparency = false, color = (:black, 0.15))

axis_real = scene_real[Axis]
axis_real[:names, :axisnames] = ("ℜ(u)", "ℑ(u)", "ℜ[D(u)]")

scene_imag = Scene(resolution = (1000, 1000))
surface!(scene_imag, vx, vy, ℑD)
xmi, ymi, zmi = minimum(scene_limits(scene_imag))
contour!(scene_imag, vx, vy, ℑD, levels = 30, linewidth = 1, transformation = (:xy, zmi))
wireframe!(scene_imag, vx_mesh, vy_mesh, ℑD_mesh, overdraw = true, transparency = false, color = (:black, 0.15))

axis_imag = scene_imag[Axis]
axis_imag[:names, :axisnames] = ("ℜ(u)", "ℑ(u)", "ℑ[D(u)]")

center!(scene_real)
scene = vbox(hbox(β_slider, w_slider, v_slider, z_slider), hbox(scene_real, scene_imag))
