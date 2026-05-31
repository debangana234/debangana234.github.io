---
layout: default
title: Blog
permalink: /blog/
---

# Blog

Intuition-first explanations of mathematics for machine learning, data science, and computational biology.

---

## Posts

<ul>
  {% for post in site.posts %}
    <li>
      <a href="{{ post.url | relative_url }}">{{ post.title }}</a>
      <br>
      <small>{{ post.date | date: "%B %d, %Y" }}</small>
    </li>
  {% endfor %}
</ul>
