---
layout: default
title: Archimedean Lattices
nav_order: 3
parent: Appendix
---
# Archimedean Lattices
The additional archimedean lattices which possess sub-hexagonal lattices include the trihexagonal, snub hexagonal, and the rhombitrihexagonal lattice (alongside their duals).
Each lattice requires a scaling factor $$\alpha_j\in\mathbb{R}$$

$$
T_j(h,k):= \alpha_j (h^2 + hk + k^2) = \alpha_j T(h,k)
$$

so that we may expression $$T_j$$ in terms of $$T(h,k)$$.

## Trihexagonal Lattice
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

$$ T_t(h,k) = \frac{\alpha_t C_T^2}{a^2} = \frac{(\frac{2}{\sqrt{3}})^2 C_T^2}{a^2} $$

so $$\alpha_t = \frac{4}{3}$$.

## Snub Hexagonal Lattice
The snub hexagonal lattice has a scaling factor $$\alpha_s = \frac{7}{3}$$. This fact can be derived by thinking of the path from the origin of our current hexagon ($$O$$), to center of the nearest hexagon ($$B$$). A path to this adjacent hexagon can be thought of in terms of the original hexagonal lattice. Moving to the center of the nearest "hexagon" $$C$$ followed by a movement of $$\frac{1}{\sqrt{3}$$ (derived in the previous lattice calculation) by an angle of $$\frac{\pi}{6}$$ we arrive at $$B$$. This gives us a triangle $$\Delta OCB$$. We know that $$OC=1$$ and $$CB=\frac{1}{\sqrt{3}}$$, and $$\angle OAC = \frac{5\pi}{6}$$, so using the law of cosines we get

$$ OB^2 = OC^2 + CB^2 - 2 OC\cdot CB = 1^2+\Big(\frac{1}{\sqrt{3}}\Big)^2-2(1)\frac{1}{\sqrt{3}}cos(\frac{\pi}{6})=\frac{7}{3} $$.

Thus $${C_T}^2=\frac{7}{3}$$ and we then get $$\alpha_s = \frac{7}{3}$$.
