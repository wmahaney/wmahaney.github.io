---
layout: default
permalink: /
redirect_from:
  - /about_me/
---

<div class="home-intro">
  {%- comment -%}
    Photo slot: add assets/img/profile.jpg (square, ~400x400) and it appears
    automatically. Until then nothing is shown.
  {%- endcomment -%}
  {%- assign photo = site.static_files | where: "path", "/assets/img/profile.jpg" | first -%}
  {%- if photo -%}
  <img class="home-photo" src="{{ photo.path | relative_url }}" alt="Photo of William E. Mahaney">
  {%- endif -%}
  <div>
    <h1 class="home-name">William E. Mahaney</h1>
    <p class="home-position">PhD student in Mathematics, Virginia Tech</p>
  </div>
</div>

<!-- TODO: replace with your 2–3 sentence research summary. -->
I work in number theory and arithmetic geometry, mainly on supersingular isogeny graphs and post-quantum cryptography.

[Email](mailto:mahaney1@vt.edu) · [GitHub](https://github.com/wmahaney) · [CV]({{ '/cv_william_mahaney.pdf' | relative_url }})

Outside of math I enjoy bouldering, Magic: The Gathering, video games, and spending time with animals.
