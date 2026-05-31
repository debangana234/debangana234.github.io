# Static GitHub Pages architecture

This version does NOT use Jekyll.

Upload the entire `Understanding-Fundamentals/` folder into the root of your `debangana234.github.io` repository.

Your repository should look like this:

debangana234.github.io/
├── index.html
└── Understanding-Fundamentals/
    ├── index.html
    ├── eigenvectors-eigendecomposition-and-pca.html
    └── assets/
        └── eigenvectors-pca/
            ├── covariance_matrix_equation.png
            ├── eigenvector_direction_visualization.gif
            ├── force_components_cos_sin.png
            ├── matrix_linear_transformation.gif
            └── pca_simple_diagram.png

Then edit the root `index.html` file and add:

<a href="/Understanding-Fundamentals/">Understanding Fundamentals</a>

or:

<a href="/Understanding-Fundamentals/">Blog</a>

After deployment, visit:

https://debangana234.github.io/Understanding-Fundamentals/

The first article will be at:

https://debangana234.github.io/Understanding-Fundamentals/eigenvectors-eigendecomposition-and-pca.html

Do not use `_posts`, `_layouts`, `blog`, or `_config.yml` for this simpler setup.
