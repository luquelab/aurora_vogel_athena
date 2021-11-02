---
layout: default
title: Archimedean Lattices
nav_order: 3
parent: Background
---
# Archimedean Lattices
The additional archimedean lattices which possess sub-hexagonal lattices include the trihexagonal, snub hexagonal, and the rhombitrihexagonal lattice (alongside their duals).
Each lattice requires a scaling factor $$\alpha_j\in\mathbb{R}$$

$$
T_j(h,k):= \alpha_j (h^2 + hk + k^2) = \alpha_j T(h,k)
$$

so that we may expression $$T_j$$ in terms of $$T(h,k)$$.

The trihexagonal lattice has a scaling factor $$\alpha_t = \frac{4}{3}$$. This comes from the fact that to navigate to an adjacent hexagon from the center of our current hexagon,
call it $$O$$, we pass through the vertex
of our current hexagon ($$B$$) rather than through the midpoint of the side ($$A$$). The distance from $$OA$$ is $$\frac{a}{2}$$ by assumption. $$OA$$ is perpendicular to $$AB$$,
which allows us to use the fact that

$$ \frac{\sqrt{3}}{2}=\sin(\frac{\pi}{3})=\frac{OA}{OB}=\frac{\frac{a}{2}}{OB}$$

which gives us

$$ OB = \frac{1}{\sqrt{3}}. $$

This means $$2 OB$$ is the distance to the nearest hexagon. The calculation of $$T$$ number boiled down to 

$$ \frac{C_T^2}{a^2} $$

thus we get

$$ T_j(h,k) = \frac{\alpha_j^2 C_T^2}{a^2} = \frac{(\frac{2}{\sqrt{3}})^2 C_T^2}{a^2} $$.
