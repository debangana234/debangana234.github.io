# Upload workflow for your existing `debangana234.github.io` site

You already have an `index.html` file. Keep it.

Upload the following new files and folders into the root of your existing `debangana234.github.io` repository:

- `_config.yml`
- `blog/`
- `_layouts/`
- `_posts/`
- `figures/`

Your final repository should look like this:

```text
debangana234.github.io/
├── index.html
├── _config.yml
├── blog/
│   └── index.md
├── _layouts/
│   ├── default.html
│   └── post.html
├── _posts/
│   └── 2026-05-31-eigenvectors-eigendecomposition-and-pca.md
└── figures/
    └── eigenvectors-pca/
        ├── covariance_matrix_equation.png
        ├── eigenvector_direction_visualization.gif
        ├── force_components_cos_sin.png
        ├── matrix_linear_transformation.gif
        └── pca_simple_diagram.png
```

Then edit your existing `index.html` and add this link somewhere visible:

```html
<a href="/blog/">Blog</a>
```

Do not create a `.nojekyll` file. Jekyll needs `_posts/` and `_layouts/` to be processed.

The blog index will be available at:

```text
https://debangana234.github.io/blog/
```

Your first post will be listed automatically on that page after GitHub Pages rebuilds.

Important:
- Do not paste the actual post into `blog/index.md`.
- Put actual blog posts inside `_posts/`.
- Post filenames must follow this pattern: `YYYY-MM-DD-title.md`.
- Equations are written with MathJax delimiters and wrapped in HTML spans/divs so Markdown does not break the math.
- GitHub's file preview may show raw math. Check the live `github.io` website URL after deployment.
