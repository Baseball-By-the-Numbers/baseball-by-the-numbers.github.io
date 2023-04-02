---
layout: post
published: true
title: Robot long jump
date: 2023/04/01
subtitle: ...
tags:
---

>**Question**: Great news! The variety of robotic competition continues to grow at breakneck pace! Most recently, head-to-head long jump contests have been all the rage.
>
>These contests consist of rounds in which each robot has a single attempt to score. In an attempt, a robot speeds down the running track (modeled as the numberline) from $0,$ the starting line, to $1,$ the takeoff point. A robot moves along this track by drawing a real number uniformly from $\left[0,1\right]$ and adding it to the robot’s current position. After each of these advances, the robot must decide whether to jump or wait. If a robot crosses the takeoff point (at $1$) before jumping its attempt receives a score of $0.$ If the robot jumps before crossing $1,$ it draws one final real number from $\left[0,1\right]$ and adds it to its current position, and this final sum is the score of the attempt.
>
>In a head-to-head contest, the two robots each have a single attempt without knowing the other’s result. In the case that they tie (typically because they both scored $0$), that round is discarded and a new round begins. As soon as one robot scores higher than the other on the same round, that robot is declared the winner!
>
>Assume both robots are programmed to optimize their probability of winning and are aware of each other’s strategies. You are just sitting down to watch a match’s very first attempt (of the first round, which may or may not end up being discarded). What is the probability that this attempt scores $0$? Give this probability as a decimal rounded to $9$ digits past the decimal point.

<!--more-->

([FiveThirtyEight](https://www.janestreet.com/puzzles/current-puzzle/))

## Solution

the jumper's decision can't depend on anything but their current position since that's all the information they have.

so, each jumper's strategy takes the form
- if current position $z$ is less than a threshold $t$, take another step
- if $z\geq t$, take the jump

each player advances down the track, step by step, until they pass $t.$ the biggest $t$ that should be considered is $1,$ or else they'll always score zero. since the steps are uniform draws from zero to $1,$ the greatest realizable score, regardless of $t,$ is $2$ (the first step can land at $z=1$, and a jump of $1$ from there would land at $2$).

## strategy outline

player $a$ wins if, eventually, they score higher than player $b.$ this can happen if

- both players don't score in the present round, and player $a$ goes on to win later, or
- player $a$ scores in the present round, but player $b$ does not.
- both players score in the present round and $s_a > s_b$,


calling the probability to not score with threshold $t$, $P_\text{zero}(t),$ and the probability of getting score $s$ given that they use threshold $t$ $P_\text{score}(s, t),$ this means

$$ 
  P(a\ \text{wins}\rvert t_a, t_b) = P_\text{zero}(t_a)P_\text{zero}(t_b)P(a\ \text{wins}\rvert t_a, t_b) + P_\text{zero}(t_b)(1-P_\text{zero}(t_a)) + 
  \int\limits_{t_b}^2 \text{d}s_b\, \int\limits_{s_b}^{2} \text{d}s_a\, P_\text{score}(s_b|t_b)P_\text{score}(s_a|t_a),
$$

or after solving for $P(a\ \text{wins})\rvert t_a, t_b)$,

$$ P(a\ \text{wins}\rvert t_a, t_b) = \dfrac{P_\text{zero}(t_b)(1-P_\text{zero}(t_a)) + 
  \int\limits_{t_b}^2 \text{d}s_b\, \int\limits_{s_b}^{2} \text{d}s_a\, P_\text{score}(s_a|t_a)P_\text{score}(s_b|t_b)}{1 - P_\text{zero}(t_a)P_\text{zero}(t_b)}.
$$

the chance to win, $P(a\ \text{wins}\rvert t_a, t_b),$ depends on both their thresholds. we can think of it as a surface above the $t_a\times t_b$ plane.

![](/img/2023-04-01-ta-tb-prob-surface.png){:width="500 px" class="image-centered"}

strategically, player $a$ should set $t_a$ so that $P(a\ \text{wins})$ has the greatest minimum with respect to $t_b.$ the same is true in the other direction. 

the game is symmetric for both players, so they'll pick the same $t=t_a=t_b$ and it suffices to find $t_a$ where

$$ 0 = \dfrac{\partial P(a\ \text{wins})}{\partial t_a}\Biggr|_{t_a=t_b}. $$

<!-- the game is symmetric for both players, so both players will pick the same $t=t_a=t_b.$ strategically, player $b$ should set $t_b$ so that $P(b\ \text{wins}|t_a, t_b)$ is maximal with respect to $t_b,$ and minimal with respect to $t_a.$ -->

to calculate the optimal threshold, we have to find expressions for $P_\text{zero}(t)$ and $P(s\rvert t).$

## finding $P_\text{zero}(t)$

if we're at position $z,$ then we can can get a zero if we step beyond $1$ on our present turn (probability $z$), or step to $z^\prime < 1$ on our present turn, and go on to get a zero from there.

$$ P_\text{zero}(z\rvert t) = z + \int\limits_z^t\text{d}z^\prime P_\text{zero}(z^\prime\rvert t). $$

taking the derivative with respect to $z,$ this becomes 

$$ P_\text{zero}^\prime(z\rvert t) = 1 - P_\text{zero}(z\rvert t), $$

which, integrating from $0$ to $z$ gets $\log\tfrac{P_\text{zero}(z\rvert t)-1}{P_\text{zero}(0\rvert t)-1} = -z. $

we want the probability of getting zero from the start of the game (when $z=0$), and we know that $P_\text{zero}(t\rvert t) = t,$ so we get

$$ P_\text{zero}(t) = 1 +(t-1) e^t. $$


## finding $P_\text{score}(s\rvert t)$

to get a non-zero score, the player first has to enter the region $\left[t, 1\right].$ the overall probability of landing inside this region is just $(1-t)$ and is uniform within it. 

say $x$ is the point where the player entered the region, and $\varepsilon$ is the size of their jump. they can jump to a point $s$ from as far away as $x=t$ or $(s-1),$ whichever is bigger. similarly, their jump can start as close as $x=1$ or $s$ itself, whichever is smaller.

with $x$ established, the size of the jump is fixed to $(s-x)$ which means that the probability of scoring $s,$ given that they score, is equal to

$$ P(s\rvert t,\ \text{score}) = \int\limits_{\max{t,(s-1)}}^{\min{s,1}} \hskip{-1em}\,dx \int\limits_0^1\,\hskip{-0.3em}d\varepsilon\, \delta(\varepsilon-(s-x)). $$

(here, $\delta()$ is the Dirac delta function)

working that out, we get

<!-- $$ 
  P(s\rvert t) = 
    \frac{1}{1-t} \begin{cases}
      (s-t) & s < 1 \\
      (1-t) & 1<s<(1+t) \\
      (2-s) & (1+t) < s.
    \end{cases}
$$ -->

$$ P(s\rvert t,\ \text{score}) = 
     \frac{1}{1-t}
\begin{cases}
    s-t & \text{if} & s < 1 \\
    1-t & \text{if} & 1<s<1+t \\
    2-s & \text{if} & 1+t < s.
\end{cases}
$$

and 

$$P_\text{score}(s\rvert t) = P(s\rvert t,\ \text{score})(1-P_\text{zero}(t)). $$

## bringing it home

putting it all together, the probability that $a$ wins, given the thresholds, is

$$ P(a\ \text{wins}\rvert t_a, t_b) = \frac{e^{t_a} (t_a-1) \left(e^{t_b} \left(-t_a^3+t_a^2-3
   (t_a-1) t_b^2+2 (t_a-1)^2 t_b+4 t_a+2
   t_b^3-6\right)+12\right)}{12 \left(\left(e^{t_a} (t_a-1)+1\right)
   e^{t_b} (t_b-1)+e^{t_a} (t_a-1)\right)}.
$$

as expected, when $t_a=t_b,$ this expression is equal to $\frac12.$

taking the derivative with respect to $t_a,$ then setting $t_a = t_b,$ we get the relatively tidy condition

$$ 
  \begin{align}
    0 &= \dfrac{\partial P(a\ \text{wins})}{\partial t_a}\Biggr|_{t_a=t_b} \\
    &= \dfrac{3 t_a-e^{t_a} (1-t_a)^2 (t_a+2)}{6 \left(e^{t_a}
   (t_a-1)+2\right) (t_a-1)}
  \end{align}
$$

with this in hand, we can root find to learn the optimal threshold, 

```mathematica
deriv = (3 ta - E^ta (-1 + ta)^2 (2 + ta))/(
  6 (2 + E^ta (-1 + ta)) (-1 + ta))
root = FindRoot[deriv, {tb, 4/10}, 
  WorkingPrecision -> 30]
  ```
  
  which gets $t_a \approx 0.416195355\ldots $

<br>
