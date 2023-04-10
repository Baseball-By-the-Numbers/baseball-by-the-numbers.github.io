---
layout: post
published: true
title: Hard choices for the Easter bunny
date: 2023/04/09
subtitle: Hoppity hop, which carton are you going to get?
tags: counting symmetry 
---

>**Question**: For Easter, you and your family decide to decorate $10$ beautiful eggs. You pull a fresh carton of eggs out of your fridge and remove $10$ eggs. There are two eggs remaining in the carton, which you return to the fridge.
>
>The next day, you open the carton again to find that the positions of the eggs have somehow changed — or so you think. Perhaps the Easter Bunny was snooping around your fridge?
>
>The $12$ slots in the carton are arranged in a six-by-two array that is symmetric upon a $180$-degree rotation, and the eggs are indistinguishable from each other. How many distinct ways are there to place two eggs in this carton? (Note: Putting two eggs in the two leftmost slots should be considered the same as putting them in the two rightmost slots, since you can switch between these arrangements with a 180-degree rotation of the carton.)
>
>**Extra credit**: Instead of two eggs remaining, suppose you have other numbers of indistinguishable eggs between zero and $12.$ How many distinct ways are there to place these eggs in the carton?

<!--more-->

([FiveThirtyEight](URL))

## Solution

if we weren't worried about rotational symmetry, then the answer would just be $\binom{12}{n},$ since that's how many ways there are to place $n$ objects in $12$ receptacles.

since we are worried about rotations, we need to make sure we count only one member of each collection rotationally equivalent arrangement.

the first thing to notice is that an arrangement with different numbers of eggs on either side of the middle can't be rotationally equivalent to itself. this makes it easy to count cases with unequal numbers of eggs. 

if there are $\ell$ eggs on the left, $r$ on the right, and $\ell \neq r,$ then there are $\binom{6}{\ell}\times\binom{6}{r}$ possible arrangements. 

### odd total number

so, for $5$ total eggs, the number of total arrangements is 

$$ \Omega(5) = \binom{6}{5}\binom{6}{0} + \binom{6}{4}\binom{6}{1} + \binom{6}{3}\binom{6}{2} = 396. $$

if we included terms with $\ell < r$ then we'd be counting the rotationally symmetric partners we want to exclude.

### even total number

when there are an even total number of eggs, it means that we have to include the case where $\ell = r.$ 

if we naively multiply the number of left and right arrangements, we will have a two-fold mess on our hands.

there are $\binom{6}{\frac12n}$ cases where an arrangement pairs with its rotationally symmetric self. these contribute one arrangement to the total. every other pairing has a rotationally symmetric duplicate represented in the product. these terms have to be downweighted by half to avoid double-counting.

this means that the proper $\ell=r$ term is

$$ \frac12\left[\binom{6}{\frac12 n}^2 - \binom{6}{\frac12 n}\right] + \binom{6}{\frac12 n} = \frac12\left[\binom{6}{\frac12 n}^2 + \binom{6}{\frac12 n}\right]. $$


### altogether

finally, for the cases beyond $n=6,$ the symmetry between holes and eggs means that $\Omega(n) = \Omega(12-n),$ and we don't have to calculate them.

$$
  \begin{array}{c|c}
     n & \Omega(n) \\ \hline
     1 & 6 \\
     2 & 36 \\
     3 & 110 \\
     4 & 255 \\
     5 & 396 \\
     6 & 472 \\
  \end{array}
$$

<br>