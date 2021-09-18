---
layout: default
title: Background
nav_order: 2
---

## Background
Put lots
of
words
in
here

### Definition
The **cross product** between two vectors $$\vec{a}$$ and $$\vec{b}$$ is defined by
<p align="center">
    $$ \vec{a} \times \vec{b} := ||\vec{a}|| ||\vec{b}||sin(\theta)\vec{n}$$
</p>
where $$\theta$$ is the angle between $$\vec{a}$$ and $$\vec{b}$$ in the plane they lie in and $$\vec{n}$$ is a unit vector perpendicular to this plane.

In 3-dimensions this is equivalent to the determinant of vector $$\vec{a}$$ and $$\vec{b}$$:
<p align="center">
    $$ \vec{a}\times\vec{b} = 
    \begin{vmatrix}
    \vec{i} & \vec{j} & \vec{k} \\
    x_1 & y_1 & z_1 \\
    x_2 & y_2 & z_2
    \end{vmatrix}$$
</p>


### Proposition
The area of the most elementary equilateral triangle contained on the hexagonal grid is denoted as $$S_0$$ and given by
<p align="center">
    $$S_0 = \frac{\sqrt{3}a^2}{4}$$
</p>
where $$a$$ is the distance between the centers of two adjacent hexagons within the hexagonal grid.


_Proof._ $$\vec{a_1} \times \vec{a_2}$$ forms a parallelogram with lengths $$\|\vec{a_1}\|$$ and $$\|\vec{a_2}\|$$. Since the basis vectors form a 60$$\degree$$ angle and both vectors have the same length, their corresponding angles must be equivalent. 

Given a triangles angles must sum to $$180\degree$$, we have $$\alpha = 60\degree$$.
Thus $$\frac{\vec{a_1} \times \vec{a_2}}{2}$$ gives the area for an equilateral triangle with side length $$a$$. 
Since
<p align="center">
    $$\|\frac{\vec{a_1} \times \vec{a_2}}{2}\| = \|\frac{\|\vec{a_1}\|\|\vec{a_2}\|sin(60\degree)\vec{n}}{2}\| = \|\frac{a^2\sqrt{3}\vec{n}}{4}\| = \frac{\sqrt{3} a^2}{4}$$
</p>
we conclude $$S_0 = \frac{\sqrt{3} a^2}{4}$$.

