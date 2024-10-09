= Matmek mandatory assignment

== 1.2.3
Want to show that 
$
u(t, x, y) = exp(iota (k_x x +  k_y y - omega t))
$
is a solution to the wave equation 
$
(diff^2 u) / (diff t^2) = c^2 nabla^2  u. 
$
We compute and find that for the above expression of u we have
$
(diff^2 u) / (diff t^2) = (- iota omega)^2 u = - omega^2 u,
$
$
(diff^2 u) / (diff x^2) = (iota k_x)^2 u = - k_x^2 u,
$
and
$
(diff^2 u) / (diff y^2) = (iota k_y)^2 u = - k_x^2 u.
$
So overall we have 
$
c^2 nabla^2 u = c^2 (- k_x^2 - k_y^2) u= - c^2 norm(k)^2 u= - omega^2 u = (diff^2 u) / (diff t^2), 
$
as we wanted to show (ie. u given by the above expression is a solution to the wave equation).

==