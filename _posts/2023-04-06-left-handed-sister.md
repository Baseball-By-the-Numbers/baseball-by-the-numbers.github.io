---
layout: post
published: false
title: 
date: 2018/04/21
subtitle:
tags:
---

>Question

<!--more-->

([FiveThirtyEight](URL))

## Solution

we want to know if the person we saw is some guy, or "the guy". the question is, how much should we increase our credence in "the guy" on the basis of hearing about his second sister.

before we hear, our relative belief is

$$
  \dfrac{P(\text{the guy}\rvert \text{sister Mary})}{P(\text{some guy}\rvert \text{sister Mary})} =   \dfrac{P(\text{sister Mary}\rvert \text{the guy})}{P(\text{sister Mary}\rvert\text{some guy})} \dfrac{P(\text{the guy})}{P(\text{some guy})}
$$
$$

and after we hear about the lefty sister, it becomes

$$
  \dfrac{P(\text{the guy}\rvert \text{sister Mary, lefty sister})}{P(\text{some guy}\rvert \text{sister Mary, lefty sister})}
$$

using bayes rule, this is equal to

$$
  \dfrac{P(\text{sister Mary, lefty sister}\rvert\text{the guy})}{P(\text{sister Mary, lefty sister}\rvert\text{some guy})}\dfrac{P(\text{the guy})}{P(\text{some guy})}
$$

using bayes rule again, this becomes

$$
  \dfrac{P(\text{lefty sister}\rvert\text{sister Mary, the guy})}{P(\text{lefty sister}\rvert\text{sister Mary, some guy})}\dfrac{P(\text{sister Mary}\rvert\text{the guy})}{P(\text{sister Mary}\rvert\text{some guy})}\dfrac{P(\text{the guy})}{P(\text{some guy})}
$$

which is just $P(\text{lefty sister}\rvert\text{sister Mary, the guy})/P(\text{lefty sister}\rvert\text{sister Mary, some guy})$ times our relative belief before hearing about the lefty sister. this update factor is called the "bayes factor" associated with the lefty sister.

the numerator is $1$ on account that "the guy" definitely has these two sisters, so we can turn our focus on the denominator $P(\text{lefty sister}\rvert\text{sister Mary, some guy}).$

<br>
