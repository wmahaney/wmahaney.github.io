---
layout: page
title: Blog
permalink: /blog/
---

<ul class="post-list">
  {%- for post in site.posts -%}
  <li>
    <span class="post-meta">{{ post.date | date: "%B %-d, %Y" }}</span>
    <h3>
      <a class="post-link" href="{{ post.url | relative_url }}">{{ post.title | escape }}</a>
    </h3>
  </li>
  {%- endfor -%}
</ul>

<p class="rss-subscribe">Subscribe <a href="{{ '/feed.xml' | relative_url }}">via RSS</a></p>
