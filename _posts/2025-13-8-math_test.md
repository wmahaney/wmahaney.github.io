---
layout: post
title: "Testing Math"
date: 2025-08-13
---

I am a mathemetician, so it is befitting I speak on it. My primary area is Number Theory, specifically supersingular elliptic curves and Post-Quantum Cryptography. I also integrate Probability Theory with Number Theory and study random properties of arithmetic data.

I also get to use this as an opportunity to test out MathJax integration in Jekyll for making my website. For the programmers, I added the block 

````html
 <!-- MathJax configuration -->
<script>
  window.MathJax = {
    tex: {
      inlineMath: [['\\(', '\\)'], ['$', '$']],
      displayMath: [['\\[', '\\]'], ['$$', '$$']]
    },
    options: {
      skipHtmlTags: ['script', 'noscript', 'style', 'textarea', 'pre', 'code']
    }
  };
</script>
````
to my _site/index.html file in my website repo. 


Now I will present a fun fact about elliptic curves and a fun fact about random variables to show off the integration. 


# Elliptic Curves Fun Fact

Let $k$ be an algebraically closed field of characteristic $p \ge 5$. An _isogeny_ of elliptic curves is a non-constant rational map $\phi: E \to \widetilde{E}$ between two elliptic curves that preserves the marked points, i.e. $\phi(o_E) = o_{\widetilde{E}}$. An _endomorphism_ of an elliptic curve $E$ is an isogeny from $E$ to itself. The _endomorphism ring_ of an elliptic curve $E$
