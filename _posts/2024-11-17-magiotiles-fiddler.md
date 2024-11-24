---
layout: post
published: true
title: My first post
date: 2024/11/23
subtitle: I'm trying this out
tags: symmetry optimization
---

Hello, world!

![](/img/table_test.png)

<!--more-->

([Fiddler on the Proof](https://thefiddler.substack.com/p/can-you-figure-out-how-magna-tiles))

## Solution

The seemingly magic behavior of the tiles tells us a lot about how magnets A-H must be oriented. Because they're impervious to rotations, magnets A, C, E, and G must have the same orientation, as must B, D, F, and H. This means we can focus on the orientations of one side. 

Imagine we stack two tiles so that A is atop A and B is atop B. This is an attractive arrangement and so must place N with S. This rules out horizontal orientations. Suppose magnet A was arranged with S on the left and N on the right. Then magnet A on the bottom would need to have the opposite arrangement to A on top, a contradiction.

![](/img/2024-11-16-tile-flip.png){:width="450 px" class="image-centered"}

### Orientation

So, the magnets must be vertically oriented. Suppose we orient magnet A so that N points up and S down. Since the interaction is unchanged when we flip one of the tiles, magnet B must have the opposite orientation to magnet A. 

![](/img/2024-11-16-tile-fields.png){:width="450 px" class="image-centered"}

Putting it all together, each tile's magnets are oriented like so:

![](/img/2024-11-17-tile-orientation-diagram.png){:width="350 px" class="image-centered"}

where U/D means the magnet's north/south end is up.

## Stable attachments

Two magnatiles snap tightly together when they're perfectly overlapped.

But if you push them out of this configuration, and slide one tile over the other, you feel attraction to a few other configurations. Likewise, if you try to slide them out of these configurations you'll feel resistance.

This is stability. A configuration is stable if the tiles resist small movements away from it. The reason a configuration is stable is that, compared to all its neighboring configurations, the attractive N-S magnet head pairs are as close together as possible while the repulsive N-N and S-S magnet head pairs are as separated as possible. 

To explore other configurations, we just need to translate and rotate one tile while keeping the other one fixed. But to find alternate configurations, we need to be able to specify them. We can specify a configuration by the displacement of the sliding tile's center of mass $\mathbf{r} = \left(x,y\right)$, and its rotation about its center of mass $\theta.$

### Setting up the tiles

Given the arrangement of magnets we figured out above, their orientations are

```mathematica
ors = { 1, -1, 1, -1, 1, -1, 1, -1 };
```

and their locations on the tile are given by

```mathematica
locs[z_]:=
     {
       {0.25, 1, z}, {0.75, 1, z}, {1, 0.75, z}, {1, 0.25, z}
     , {0,75, 0, z}, {0.25, 0, z}, {0, 0.25, z}, {0, 0.75, z}
     };
```

where $z$ keeps track of the vertical separation between the planes of the two tiles.

### Rotations and translations

As we said above, sliding the tile around amounts to spinning and translating all the coordinates of the tile:

```mathematica
locsRotAndTrans[x_, y_, theta_, locs_]:=(
   COM = {0.5, 0.5, 0};

   locsCOM = Table[
               locs[[j]] - COM
               , {j, 1, 8}];

   locsCOMRot = Table[
                  RotationMatrix[theta, {0, 0, 1}] . locsCOM[[j]]
                  , {j, 1, 8}];

   locsRot = Table[
               locsCOMRot[[j]] + COM
               , {j, 1, 8}];

   locsRotTrans = Table[
                     locsRot[[j]] + {x, y, 0}
                     , {j, 1, 8}];

   Return[locsRotTrans];
)
```

The interaction between two tiles is the sum of the interactions between the component magnets on either tile. 

$$ E = \sum\limits_{i,j\in\{A\,\ldots,H\}} \frac{o_io_j}{\ell_{ij}} $$

where the $o_i$ are the orientations we defined above, and $\ell_{ij}$ is the displacement between the tile locations we defined above:

```mathematica
   locsRT = rotAndTransLocs[x, y, theta, locs[-0.05]];
   interaction =
      Sum[
         ors[[i]] * ors[[j]] / Norm[locs[0.05][[i]] - locsRT[[j]]]
         , {i, 1, 8}
         , {j, 1, 8}
      ];
```

### Exploring configuration space

Doing all that gets us the two tile's interaction energy $E$ parameterized in terms of $x$, $y$, and $\theta$: $E(x,y,\theta).$ This forms a landscape of energies where depressions in the landscape correspond to the stable arrangements we're trying to find. 

For example, if we hold $\theta$ fixed at zero and vary the position of the second tile's center of mass, the energy landscape looks like

![](/img/2024-11-18-surface-blue-a.png){:width="450 px" class="image-centered"}

The deep well in the middle is the configuration where the tiles overlap exactly, and it's surrounded by shallower wells that correspond to lesser stable states.

If, instead, we hold $\theta$ fixed at $\pi/4$ and vary the center of mass, we see a repulsive peak in the center that's surrounded by shallow wells 

![](/img/2024-11-18-surface-blue-b.png){:width="450 px" class="image-centered"}

However, these are snapshots at specific values of $\theta$, mere projections of the full space. In general, we can find the stable configurations by following the gradient downhill — start at a random initial position $(x_0, y_0, \theta_0)$ and then take small steps in the direction of steepest descent:

$$
   \langle x, y, \theta\rangle_t = \langle x, y, \theta\rangle_{t-1} - \eta \nabla E(x, y, \theta).
$$

Since there are multiple minima, we'll have to run gradient descent many times, from random starting positions, until we're pretty sure we've found them all (the common ones at least).

### Behold the configurations

Carrying this out, we find the four confgurations below:

![](/img/2024-11-17-grid-no-offset.png){:width="450 px" class="image-centered"}

The first is the ultra stable configuration with energy $-54.27$. The second features the red tile pushing a side across the corner of the purple tile until it's just past the purple diagonal. This places two pairs of opposite monopoles in close proximity for an overall energy of $-13.31$. This is followed closely in third place by an arrangment where adjacent sides overlap so that two monopoles are also in close proximity with overall energy $-12.74$. Finally, there's an arrangement that's essentially the second had we stopped sliding the red tile early. This configuration places two opposite monopoles in close proximity, but loses the stabilizing long distance attractive interactions of the first two, yielding an overall energy of $-9.34$.

As it happens, these closely match the stable arrangements you find by sliding two magnatiles in real life. It is really neat how well the model works given the monopole approximation.

### Meet the real number two

With real magnatiles you'll also find that there's an even more stable arrangement formed by placing two magnatiles side by side. However, as constructed, our model cannot find or quantify it. The first reason is that our model stipulates the $z$-separation of the tile to be nonzero. The second is that our model puts the monopoles directly at the edge which means there'd be zero distance between them when the tiles are placed side by side. To accomodate this configuration, we made a modified model where the monopoles are pushed in by a small amount ($\delta = 0.03$) from the edge. 

When this is done, we find the original stable arrangements with approximately the same energies:

![](/img/2024-11-17-grid-offset.png){:width="450 px" class="image-centered"}

but we can also quantify the stability of the side by side arrangement, which jumps into second place with $E = -25.57$:

![](/img/2024-11-17-grid-offset-side.png){:width="300 px" class="image-centered"}

<br>
 
