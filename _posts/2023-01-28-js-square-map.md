---
layout: post
published: true
title: Chain of squares 
date: 2023/01/28
subtitle: How long can you keep your square away from the covetous origin?
tags: fixed-point dynamics 
---

>**Question**:
>Assign four nonnegative integers to the corners of a square, which we designate the active square. During a step, for each side of the active square, the absolute difference between the numbers on that side’s endpoints is assigned to its midpoint. Then these four new midpoints are connected into a new square (tilted $45$ degrees from the previous). This new smaller square becomes the active square. Continue these steps until the active square has all zeroes on its corners.
>
>Define $f(a, b, c, d)$ to be the total number of squares drawn during this process when beginning with the numbers $(a, b, c, d)$ written on the starting square in clockwise order. For example, given a starting arrangement of $(10, 6, 3, 1),$ we would get the sequence of
>
>$(4, 3, 2, 9)$
>$(1, 1, 7, 5)$
>$(0, 6, 2, 4)$
>$(6, 4, 2, 4)$
>$(2, 2, 2, 2)$
>$(0, 0, 0, 0)$
>
> where the game ends (pictured above). So $f(10, 6, 3, 1) = 7.$ And trivially, $f(0, 0, 0, 0) = 1.$
>
> Consider the set $S = \\{(a, b, c, d) \rvert a, b, c, d\in\mathbb{Z}, 0 \leq a, b, c, d \leq 10^7\\}.$ 
>
> Let $M$ be the maximum value $f$ obtains on $S.$ Find $(a, b, c, d)$ in $S$ with minimum sum $(a+b+c+d)$ where $f(a, b, c, d) = M.$

<!--more-->

([Jane Street](https://www.janestreet.com/puzzles/lesses-more-index/))

## Solution

first of all, the map from one square to the next is 

$$ \left(a,b,c,d\right) \overset{M}{\longrightarrow} \left(\lvert a-b\rvert, \lvert b-c\rvert, \lvert c-d\rvert, \lvert d-a\rvert\right). $$

<!-- from this we can see two properties of the map:

1. the map is  -->

we can follow the map on a couple of squares to get the hang of it

![drawing of squares](/img/2023-01-30-square-drawing.png){:width="450 px" class="image-centered"}

we want to find the starting configuration $(a,b,c,d)$ that takes the most map iterations to reach the $\left(0,0,0,0\right)$ square. 

from the map, we can see two properties of $f$:

1. $f$ is invariant under uniform translations

$$ f(a,b,c,d) = f(a+k,b+k,c+k,d+k) $$
  
2. $f$ is invariant to uniform scaling

$$ f(a,b,c,d) = f(\gamma a, \gamma b, \gamma c, \gamma d). $$

the second property shows that we can think of $(a,b,c,d)$ as a relative composition (normalize it so that its components sum to $1$). in principle, if we can find a composition $(a,b,c,d)$ that maps back to itself, the value of $f$ would be infinite. 

we will look for such a composition, and then look for the closest integer approximation we can find. the bigger numbers we use, the more closely we should expect to get.

if $(a,b,c,d)$ maps to itself, then, up to an overall multiplicative constant:

$$ (a,b,c,d) = \left(\lvert a-b\rvert, \lvert b-c\rvert, \lvert c-d\rvert, \lvert d-a\rvert\right). $$

since the two compositions are the same, ratios between like components are the same. taking the ratio of the first and second components of each composition, we get:

$$
  \begin{align} 
      \frac{a}{b} &= \dfrac{\lvert a - b\rvert}{\lvert b-c\rvert} \\
    &= \frac{a}{b}\dfrac{\lvert 1-\frac{b}{a}\rvert}{\lvert 1-\frac{c}{b}\rvert}.
  \end{align}
$$

which shows that $b/a = c/b.$ the same comparison for the second and third components shows that $c/b = d/c$ as well.

putting it together, we have shown that $a = \gamma b = \gamma^2 c = \gamma^3 d,$ which shows that $(a,b,c,d)$ is monotonic

now, if we compare the first and fourth components, we get (assuming $a>b>c>d,$ and therefore $\gamma > 1$):

$$
  \begin{align}
    \frac{d}{a} &= \dfrac{\lvert d-a\rvert}{\lvert a-b\rvert} \\
    \gamma^3 &= \dfrac{1 - \gamma^3}{1-\gamma},
  \end{align}
$$

or $2\gamma^3 -\gamma^4 - 1 = 0$ which has one relevant root

$$
  \begin{align}
    \gamma &= \frac13 \left(1 + \sqrt[3]{\left(19 - 3 \sqrt{33}\right)} + \sqrt[3]{19 + 3 \sqrt{33}}\right) \\
    &\approx 1.8393
  \end{align}
$$

so, an ideal composition vector is given by 

$$ \dfrac{1}{1+\gamma+\gamma^2+\gamma^3}\left(1,\gamma,\gamma^2,\gamma^3\right) $$ 

which is stationary under the map, as expected. from here we can just scaled the composition vector, looking for its closest integer approximation, which is

$$ \vec{\phi} = \left(1655616,3045153,5600910,10301680\right). $$

running this through $f$ gets $f(\vec{\phi}) = 44.$

because of the first property, we can subtract the minimum entry from each component of the composition vector without changing $f,$ so the vector with minimum sum is 

$$ \vec{\phi} = \left(0,1389537,3945294,8646064\right). $$

<br>
