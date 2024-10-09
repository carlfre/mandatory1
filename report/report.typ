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

== 1.2.4
We have the following discretization
$
(u^(n+1)_(i, j) - 2 u^n_(i,j) + u^(n-1)_(i,j)) / (Delta t^2)
 = c^2 [ (u^n_(i+1, j) - 2u^n_(i,j) + u^n_(i-1, j)) / h^2 + (u^n_(i, j+1) - 2u^n_(i, j) + u^n_(i, j-1) ) / (h^2) ]
$
Rearranging, we find:
$
(u^(n+1)_(i, j) - 2 u^n_(i,j) + u^(n-1)_(i,j)) &= C^2 [ (u^n_(i+1, j) - 2u^n_(i,j) + u^n_(i-1, j))  + (u^n_(i, j+1) - 2u^n_(i, j) + u^n_(i, j-1) )  ] \
&= C^2 [ u^n_(i+1, j)  + u^n_(i-1, j)  + u^n_(i, j+1)  + u^n_(i, j-1) -  4 u^n_(i, j)  ]
$
We insert for $
u^n_(i,j) = exp(iota (k h (i + j) - tilde(omega) n Delta t))
$
and divide by $u_(i,j)^n$ and get 

$
exp(-iota tilde(omega) Delta t) - 2 + exp(iota tilde(omega) Delta t) = C^2 [ 2 exp(iota k h) + 2 exp(- iota k h) - 4].
$
This further rewrites to:
$
2cos(tilde(omega) Delta t ) - 2 = C^2 [4 cos(k h) - 4 ].
$
We now assume $C = 1 / sqrt(2)$. dividing both sides by 2, we get:
$
cos(tilde(omega) Delta t ) - 1 = cos(k h) - 1,
$
so overall we have the condition 
$
cos(tilde(omega) Delta t) = cos(k h).
$
we have a solution for 
$
tilde(omega) Delta t = k h, 
$
ie.
$
tilde(omega) = (k h) / (Delta t) = k c / C.
$
Furthermore, we have 
$
norm(arrow(k)) = sqrt(k_x^2 + k_y^2) = sqrt(2) k.
$
Using this and again inserting for $ C = 1 / sqrt(2)$, we find that 
$
tilde(omega) = norm(arrow(k)) c = omega. 
$

== 1.2.6
Check the readme in the top folder!