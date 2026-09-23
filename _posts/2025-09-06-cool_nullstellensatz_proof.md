---
layout: post
title: "A Cool Nullstellensatz Proof"
date: 2025-09-06
---

Today I was sitting in my Algebraic Geometry class and we were discussing Hilbert's Nullstellensatz. I have seen proofs of it before using Noether normalization, but today we saw a slick proof of it I thought was particularly satisfying. I wanted to share it here. I thank Leonardo Constantin Mihalcea for the proof.

The one piece of input we take as given is **Zariski's lemma**: if a field $$K$$ is finitely generated as an algebra over a field $$k$$, then $$K$$ is a finite (hence algebraic) extension of $$k$$. This is where the real depth of the Nullstellensatz lives; the argument below shows how quickly everything else follows from it.

## Weak Nullstellensatz

Let $$k$$ be an algebraically closed field. Then

1) All maximal ideals of $$k[x_1, \ldots, x_n]$$ are of the form $$(x_1 - a_1, \ldots, x_n - a_n)$$ for some $$a_i \in k$$.

2) For all proper ideals $$\f{a} \subset k[x_1, \ldots, x_n]$$, the variety $$V(\f{a})$$ is non-empty.

Proof:

1) If $$\f{m}$$ is a maximal ideal then $$k[x_1, \ldots, x_n]/\f{m}$$ is a field that is finitely generated as a $$k$$-algebra, so by Zariski's lemma it is an algebraic field extension of $$k$$. As $$k$$ is algebraically closed this field extension is trivial, so $$k[x_1, \ldots, x_n]/\f{m} \cong k$$. Then $$\f{m}$$ is the kernel of the surjective map $$k[x_1, \ldots, x_n] \to k$$ sending $$x_i \mapsto a_i$$ for some $$a_i \in k$$. Thus $$\f{m} = (x_1 - a_1, \ldots, x_n - a_n)$$.

2) If $$\f{a}$$ is a proper ideal then $$\f{a}$$ is contained in some maximal ideal $$\f{m}$$. By part 1 we have $$\f{m} = (x_1 - a_1, \ldots, x_n - a_n)$$ for some $$a_i \in k$$. By the order reversing nature of the $$V(\cdot)$$ map we have

$$
\{ (a_1, \ldots, a_n)\} = V(\f{m}) \subset V(\f{a}),
$$

which gives the desired conclusion. $$\square$$

## Strong Nullstellensatz

Thm: (Hilbert's Nullstellensatz) Let $$k$$ be an algebraically closed field and $$\f{a} \subset k[x_1, \ldots, x_n]$$ be an ideal. Then $$I(V(\f{a})) = \sqrt{\f{a}}$$.

Proof: If $$\f{a} = k[x_1, \ldots, x_n]$$ then $$V(\f{a}) = \emptyset$$ and $$I(V(\f{a})) = k[x_1, \ldots, x_n] = \sqrt{\f{a}}$$. So assume $$\f{a}$$ is a proper ideal.

We know already that $$\sqrt{\f{a}} \subset I(V(\f{a}))$$, so we need to show the reverse inclusion. Let $$f \in I(V(\f{a}))$$. If $$f = 0$$ then $$f \in \sqrt{\f{a}}$$ trivially, so assume $$f \neq 0$$.

By Hilbert's Basis Theorem $$\f{a}$$ is finitely generated, say $$\f{a} = (f_1, \ldots, f_r)$$. The slick step, known as the **Rabinowitsch trick**, is to add a new variable $$y$$ and consider the ideal $$\f{b} = (f_1, \ldots, f_r, 1 - yf) \subset k[x_1, \ldots, x_n, y]$$. If $$\f{b}$$ were a proper ideal then by the weak Nullstellensatz there would be some $$(a_1, \ldots, a_n, b) \in k^{n+1}$$ such that $$f_i(a_1, \ldots, a_n) = 0$$ for all $$i$$ and $$1 - b f(a_1, \ldots, a_n) = 0$$. But $$(a_1, \ldots, a_n) \in V(\f{a})$$ and $$f \in I(V(\f{a}))$$, so $$f(a_1, \ldots, a_n) = 0$$ and hence $$1 - b \cdot 0 = 1 \neq 0$$, a contradiction. So $$\f{b}$$ must be the unit ideal. Thus there are polynomials $$g_1, \ldots, g_r, h \in k[x_1, \ldots, x_n, y]$$ such that

$$
1=f_1 g_1 + \cdots + f_r g_r + (1 - yf)h.
$$

Applying the ring homomorphism $$k[x_1, \ldots, x_n, y] \to k(x_1, \ldots, x_n)$$ that fixes each $$x_i$$ and sends $$y \mapsto 1/f$$ gives

$$
1 = f_1 g_1(x_1, \ldots, x_n, 1/f) + \cdots + f_r g_r(x_1, \ldots, x_n, 1/f),
$$

since the last term becomes $$(1 - f/f)\,h(x_1, \ldots, x_n, 1/f) = 0$$. Multiplying by a sufficiently high power $$f^m$$ to clear denominators gives $$f^m = f_1 h_1 + \cdots + f_r h_r$$ for some $$h_i \in k[x_1, \ldots, x_n]$$. Thus $$f^m \in \f{a}$$, so $$f \in \sqrt{\f{a}}$$. $$\square$$
