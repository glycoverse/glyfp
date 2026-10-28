# Logo source and reproduction

`man/figures/logo.svg` is the editable source, created for glyfp. The PNG is a
transparent export for README and package documentation.

References visually inspected: glyenzy (cream hexagon, icon above wordmark),
glyrepr (black hexagon, structure and bits), glymotif (green hexagon, glycan search),
and glyhub's editable SVG (424 x 489 geometry, 12-unit border, 50-unit wordmark).
No separate logo style guide was found in the inspected package directories.

The teal palette (#416B69 / #92C7C2) distinguishes glyfp while retaining the
Glycoverse hexagon, flat icon, dark lower wordmark, and italic function suffix.
White fingerprint ridges surround a stylized SNFG-inspired branch; this is a
brand symbol, not a chemically specified example glycan. Artwork is original,
with no external raster, font, or model assets embedded.

Rebuild from the package root with librsvg:

```sh
rsvg-convert -o man/figures/logo.png man/figures/logo.svg
```

Keep the SVG as the authoritative source. Both assets use the package MIT license.
