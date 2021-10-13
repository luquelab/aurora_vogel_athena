---
layout: default
title: Background
nav_order: 3
has_children: true
---

## Background

### Definition
>
> $$ \vec{a} \times \vec{b} := ||\vec{a}|| ||\vec{b}||sin(\theta)\vec{n}$$
>
> where $$\theta$$ is the angle between $$\vec{a}$$ and $$\vec{b}$$ in the plane they lie in and $$\vec{n}$$ is a unit vector perpendicular to this plane.
>
> In 3-dimensions this is equivalent to the determinant of vector $$\vec{a}$$ and $$\vec{b}$$:
>
> $$ \vec{a}\times\vec{b} = 
> \begin{vmatrix}
> \vec{i} & \vec{j} & \vec{k} \\
> x_1 & y_1 & z_1 \\
> x_2 & y_2 & z_2
> \end{vmatrix}$$


### Proposition
> The area of the most elementary equilateral triangle contained on the hexagonal grid is denoted as $$S_0$$ and given by
>
> $$S_0 = \frac{\sqrt{3}a^2}{4}$$
>
> where $$a$$ is the distance between the centers of two adjacent hexagons within the hexagonal grid.


_Proof._ $$\vec{a_1} \times \vec{a_2}$$ forms a parallelogram with lengths $$\|\vec{a_1}\|$$ and $$\|\vec{a_2}\|$$. Since the basis vectors form a 60$$^\circ$$ angle and both vectors have the same length, their corresponding angles must be equivalent. 

Given a triangles angles must sum to $$180^\circ$$, we have $$\alpha = 60^\circ$$.
Thus $$\frac{\vec{a_1} \times \vec{a_2}}{2}$$ gives the area for an equilateral triangle with side length $$a$$. 
Since

$$\|\frac{\vec{a_1} \times \vec{a_2}}{2}\| = \|\frac{\|\vec{a_1}\|\|\vec{a_2}\|sin(60^\circ)\vec{n}}{2}\| = \|\frac{a^2\sqrt{3}\vec{n}}{4}\| = \frac{\sqrt{3} a^2}{4}$$

we conclude $$S_0 = \frac{\sqrt{3} a^2}{4}$$. 
$$\hspace{720px} \blacksquare$$


We can calculate the area of any equilateral triangle contained in the hexagonal grid as the number of base triangles contained. We call this the $$T$$-number of such a triangle.


### Theorem
>The $$T$$-number of any triangle is given by the formula
>$T = h^2 + hk + k^2$


_Proof._ Assume some vector $$\vec{C_T}=(h,k)$$. Denote it's length by $$\|C_T\|$$. This will be the length of the equilateral triangle we consider.
Notice that

$$T = \|\frac{\vec{C_T}\times\vec{C_T}}{2 }\|\frac{1}{S_0} = \frac{\frac{\sqrt{3}}{4} c_T^2}{\frac{\sqrt{3}}{4}a^2} = \frac{c_T^2}{a^2}$$

By the law of cosines we have


$$\begin{align*}
    c_T^2 &= (ha)^2 + (ka)^2 - 2(ha)(ka)cos(120^\circ) \\
    &= a^2(h^2 + k^2 - 2hk(-\frac{1}{2})) \\
    &= a^2(h^2 + k^2 + hk)
\end{align*}$$


Thus $$T = \frac{c_T^2}{a^2} = h^2 + hk + k^2$$. 
$$\hspace{720px} \blacksquare$$


**_Remark_** Although $$S_0$$ divides every equilateral triangle evenly, it is not necessarily the case that every equilateral triangle can be evenly broken into these base triangles.


Each triangle's $$T$$-number is solely determined by the length of $$\vec{C_T}$$ since each triangle is equilateral. Thus it is no surprise that a $$T$$-number is not unique to a particular $$(h,k)$$. An example is $$(7,0)$$ and $$(5,3)$$ which share $$T=49$$. In order to classify triangles like these which share certain geometrical properties (base triangle size),  Caspar and Klug proposed a reorganization of the $$T$$-structures in terms of classes $$P$$ defined as follows:
### Definition
A $$P$$ number of a $$T$$-structure is defined as

$$\begin{align*}
    P = h_0^2+h_0 k_0 + k_0^2
\end{align*}$$

where $$h = h_0 f$$, $$k = k_0 f$$ and $$f = gcd(h,k)$$.


## Elongated Structures
There are three different elongated structures which we can construct: two-fold, three-fold, and five-folds. We can use the symmetry number of such structures to conclude how many different types of triangles must be present within the body. This is calculated by multiplying the axial symmetry by the symmetry of the caps. The 5-fold case has a symmetry number of $$5*2 = 10$$.

$$\begin{array}{||c c c c ||}
 \hline
 \text{Axial symmetry} & \text{Symmetry number} & \text{Num of body} \Delta & \text{Types of body} \Delta \\
 \hline
 \hline
 5\text{-fold} & 10 & 10 & 1 \\ 
 \hline
 3\text{-fold} & 6 & 12 & 2 \\
 \hline
 2\text{-fold} & 4 & 12 & 3 \\
 \hline
\end{array}$$

It will be convenient to work with another pair of coordinates on our hexagonal lattice. Let $$\vec{a_1}'$$ and $$\vec{a_2}'$$ be vectors $$\vec{a_1}$$ and $$\vec{a_2}$$ rotated counter-clockwise by $$60^\circ$$.

We need a way to convert hexagonal coordinates to cartesian coordinates in order to interpret the cross product of the vectors in hexagonal coordinates. Thus we have the following formulas.

### Proposition
> $$(h\vec{a_1},k\vec{a_2})_H=(h+k\cos(\frac{\pi}{3}),k\sin (\frac{\pi}{3}))_C$$ 
> and similarly 
> $$(h'\vec{a_1}',k'\vec{a_2}')_H=(h'\cos(\frac{2\pi}{3})+k'cos(\frac{2\pi}{3}),h'\sin(\frac{2\pi}{3})+k'\sin(\frac{2\pi}{3}))_C$$.

_Proof._ Consider $$(h\vec{a_1},k\vec{a_2})_H$$. Let us treat the $$X$$-axis as aligned with the $$\vec{a_1}$$ axis. Thus an $$h$$ step in the $$\vec{a_1}$$ direction gives an equivalent value: $$ha$$. A $$k$$ step in the $$\vec{a_2}$$ direction forms a right triangle with line created by the projection of $$\vec{a_2}$$ onto the $$X$$-axis. Thus we have
$$
\begin{align*}
\cos(\frac{\pi}{3})=\frac{x}{k} \implies x=\frac{k}{2}
\sin(\frac{\pi}{3})=\frac{y}{k} \implies y=\frac{\sqrt{3}k}{2}
\end{align*}
$$
where $$x,y$$ correspond to the legs of the right triangle which align with the $$X$$ and $$Y$$ axes. Thus we have in total that
$$
\begin{align*}
(h\vec{a_1},k\vec{a_2})_H=(h+k\cos(\frac{\pi}{3}),k\sin (\frac{\pi}{3}))_C
\end{align*}
$$


### Definition
> Define the **body-vector $$\vec{C_Q}$$** as the vector which connects a pentamer from one cap to the closest pentamer on the opposing cap
> $$\vec{C_Q} := h' \vec{a_1}' + k' \vec{a_2}'$$
We can then define $$T_1$$ as the area of the body triangle (made up by sides $$\vec{C_Q}$$ and $$\vec{C_T}$$ divided by the area of the basic triangles (so it is analogous to the $$T$$-number but for the body triangle). We can now prove the following formula for $$T_1$$.

_Proof._ Assume vector $$\vec{C_T}=(h,k)$$ and $$\vec{C_Q}=(h',k')$$. 
