---
layout: post
published: false
title: Multiverse Carpentry
date: 2021/02/28
---

>**Question**: a regular staircase is built from blocks and the blocks in each level are different colors. The staircase can be built in whatever order that's physically possible. In addition, the surface the stairs are built on is slightly sloped, so that any placed block slides forward until it hits the wall. How many possible ways are there to build a $4$-level staircase? An $n$-level staircase?

<!--more-->

([FiveThirtyEight](URL))

## Solution

### Getting started

I got started by building some staircases. 

For the two level staircase, there's no choice for the first block, we **have** to start with an "A" block, since every other block needs a foundation on which to be placed. Once the first "A" block is placed, we have a choice for what to do second, we can either place the "B" block on top of the first "A" block, or we can place the second "A" block and then place the "B" block last. This makes for a total of $2$ ways to build the $N = 2$ level staircase.

Moving to the three level staircase, we again need to kick things off with an "A" block. From there we can either place a "B" block (on top of the first "A" block) or we can place another "A" block, extending the first level. In the situation where we placed a "B" block second, we now have the option to place a "C" block, or to place another "A" block. 

### Some observations

You can carry on like this, but it gets out of control fast. At this point I noticed a couple of patterns. 

1. we can always place an unplaced "A" block.
2. we can place a "B" block if the number of unplaced "A" blocks is equal to or less than the number of placed "B" blocks (and the same condition for "C" compared to "B").
3. if multiple blocks have legal placements, then form a branch for each possibility.
4. if we have filled out column, say by placing "A-B-C" in a stack, then the remaining $2$ "A" blocks and $1$ "B" block present the same problem as if we had the $2$-stair problem. 
5. if we have $1$ block left, then there's only one way to place it.

On account of the Observation $3$, I started using a notation in hopes to find a recursion: $\Omega(a,b,c)$ is the number of ways that $a$ remaining "A" blocks, $b$remaining  "B" blocks, and $c$ remaining "C" blocks can be placed in the staircase.

### Using the notation

Let's use the notation to recreate the two stairs result, and get a number for the three stairs result.

**Two stairs**

Starting with $2$ "A" blocks, and $1$ "B" block, we are trying to calculate $\Omega(2,1).$ With no blocks placed, we are obligated to place an "A", so we get

$$ \Omega(2,1) = \Omega(1,1).$$

Since $a = b$ we are free to place an "A" or a "B". According to Observation $2$, this means that

$$ \Omega(1,1) = \Omega(1,0) + \Omega(0,1). $$

According to Observation $5$, $\Omega(1,0)$ and $\Omega(0,1)$ are both equal to $1.$ Altogether, this shows that $\Omega(2,1) = 2.$

**Three stairs**

Now let's do the three stair case, moving a little more quickly this time. We start with $\Omega(3,2,1)$



<br>