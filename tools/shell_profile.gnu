#!/usr/bin/env gnuplot
# Smoothly truncated power-law shell profile: implementation of Sec. 7 of
# doc/shell_profile.ipynb.
#
#   rho(r) = 0 | parabola | r^alpha | parabola | 0,   C^1 at both joints
#
# Usage:  ./shell_profile.gnu     (writes shell_profile.png)

# ---- control parameters ----------------------------------------------------
# Radii are in units of the outer radius of the shell (the shell ends at r=1).
#
#   h      inner hole radius: the density is zero for r < h         (0 <= h)
#   q      rim width: the density rises from zero to the power law
#          on [h, h+q] and falls back to zero on [1-q, 1]     (h+q < 1-q)
#   alpha  power-law index of the body of the shell, rho ~ r^alpha,
#          on [h+q, 1-q]; the profile is non-negative as long as
#          -2*(1-q) <= alpha*q <= 2*(h+q)
#
h     = 0.1681
q     = 0.1
alpha = -6.0
# -----------------------------------------------------------------------------

a = h + q          # joints of the rims with the power law
b = 1.0 - q
if (h < 0 || a >= b) { print "ERROR: need 0 <= h < h+q < 1-q"; exit }

fourpi = 4.0*pi
f(r)  = r**alpha                 # power law
df(r) = alpha*r**(alpha - 1.0)   # and its derivative

# inner rim: P- = am*(r-h)^2 + bm*(r-h)
am = (q*df(a) - f(a))/q**2
bm = (2.0*f(a) - q*df(a))/q
# outer rim: P+ = ap*(r-1)^2 + bp*(r-1)
ap = -(f(b) + q*df(b))/q**2
bp = -(2.0*f(b) + q*df(b))/q

# antiderivatives of r^2 * P-(r) and r^2 * P+(r), with s = r-h and t = r-1
F(s) = am/5*s**5 + (bm + 2*h*am)/4*s**4 + (2*h*bm + h*h*am)/3*s**3 + h*h*bm/2*s**2
G(t) = ap/5*t**5 + (bp + 2*ap)/4*t**4 + (2*bp + ap)/3*t**3 + bp/2*t**2

# mass of the power-law body between a and r (alpha = -3 is a special case)
is_log = (abs(alpha + 3.0) < 1e-12)
body(r) = is_log ? fourpi*log(r/a) : fourpi*(r**(alpha+3) - a**(alpha+3))/(alpha+3)

# integration constants: mass inside the joints, and the total mass
M_a   = fourpi*F(q)
M_b   = M_a + body(b)
M_tot = M_b - fourpi*G(-q)
rho0  = 1.0/M_tot     # normalize to unit total mass

# density, its derivative and the enclosed mass
rho(r)  = rho0*(r < h ? 0 : r < a ? am*(r-h)**2 + bm*(r-h) \
              : r < b ? f(r) : r <= 1 ? ap*(r-1)**2 + bp*(r-1) : 0)
drho(r) = rho0*(r < h ? 0 : r < a ? 2*am*(r-h) + bm \
              : r < b ? df(r) : r <= 1 ? 2*ap*(r-1) + bp : 0)
mass(r) = rho0*(r < h ? 0 : r < a ? fourpi*F(r-h) \
              : r < b ? M_a + body(r) \
              : r <= 1 ? M_b + fourpi*(G(r-1) - G(-q)) : M_tot)

# non-negativity criteria of Sec. 5
km = alpha*q/a
kp = alpha*q/b
if (km > 2.0 || kp < -2.0) {
  print sprintf("WARNING: the profile is negative inside a rim (kappa- = %g, kappa+ = %g)", km, kp)
}
print sprintf("h = %g, q = %g, alpha = %g: rho0 = %g, M(h+q) = %g, M(1-q) = %g, M_tot = %g", \
              h, q, alpha, rho0, mass(a), mass(b), mass(1.0))

# the bare power law the rims are tangent to, continued through the rims
powerlaw(r) = (r >= h && r <= 1) ? rho0*f(r) : NaN

# ---- plot ------------------------------------------------------------------
set terminal pngcairo size 800,560 font ",12"
set output "shell_profile.png"

set samples 4000
set xrange [0:1.15]
# maximum of the density: at one of the joints, or at the vertex of a rim
# parabola if it lies inside the rim
max2(x, y) = x > y ? x : y
sm = (am < 0) ? -bm/(2*am) : -1     # vertex of P-, s = r-h
tp = (ap < 0) ? -bp/(2*ap) :  1     # vertex of P+, t = r-1
rho_max = max2(rho(a), rho(b))
if (sm > 0 && sm < q)  { rho_max = max2(rho_max, rho(h + sm)) }
if (tp < 0 && tp > -q) { rho_max = max2(rho_max, rho(1 + tp)) }
set yrange [1e-4*rho_max:1.1*rho_max]
set xlabel "r"
set ylabel "{/Symbol r}(r)"
set title sprintf("shell profile: h = %g, q = %g, {/Symbol a} = %g, unit total mass", h, q, alpha)
set grid lc rgb "#dddddd"
set key inside top center
set logscale y

# rims (shaded) and their edges
set style rect fc rgb "#d9d9d9" fs solid noborder behind
set obj 1 rect from h, graph 0 to a, graph 1
set obj 2 rect from b, graph 0 to 1, graph 1
set arrow 1 from h, graph 0 to h, graph 1 nohead dt 2 lc rgb "#999999"
set arrow 2 from a, graph 0 to a, graph 1 nohead dt 2 lc rgb "#999999"
set arrow 3 from b, graph 0 to b, graph 1 nohead dt 2 lc rgb "#999999"
set arrow 4 from 1, graph 0 to 1, graph 1 nohead dt 2 lc rgb "#999999"
set label 1 "rim" at h + q/2, graph 0.96 center tc rgb "#666666"
set label 2 "rim" at 1 - q/2, graph 0.96 center tc rgb "#666666"

plot rho(x)      w l lw 2 lc rgb "#3b528b" t "{/Symbol r}(r)", \
     powerlaw(x) w l lw 1 dt 3 lc rgb "#3b528b" t "power law r^{/Symbol a}"
