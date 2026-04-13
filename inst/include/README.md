# inst/include — Vendored C++ Header Libraries

This directory contains header-only C++ libraries vendored directly into the
package so no system installation is required.

## Contents

| Directory    | Library   | Version | Source                                      |
|--------------|-----------|---------|---------------------------------------------|
| `xtl/`       | xtl       | 0.7.7   | https://github.com/xtensor-stack/xtl        |
| `xtensor/`   | xtensor   | 0.25.0  | https://github.com/xtensor-stack/xtensor    |

## Populating the headers

Run the vendoring script from the repository root:

```bash
bash tools/vendor-headers.sh
```

This downloads the correct tagged releases of xtl and xtensor and extracts
only the header files into this directory.

## Manual vendoring

### xtl 0.7.7

```bash
curl -L https://github.com/xtensor-stack/xtl/archive/refs/tags/0.7.7.tar.gz \
  | tar -xz --strip-components=2 -C inst/include/xtl \
       'xtl-0.7.7/include/xtl'
```

### xtensor 0.25.0

```bash
curl -L https://github.com/xtensor-stack/xtensor/archive/refs/tags/0.25.0.tar.gz \
  | tar -xz --strip-components=2 -C inst/include/xtensor \
       'xtensor-0.25.0/include/xtensor'
```

## License

- **xtl**: BSD 3-Clause (see https://github.com/xtensor-stack/xtl/blob/master/LICENSE)
- **xtensor**: BSD 3-Clause (see https://github.com/xtensor-stack/xtensor/blob/master/LICENSE)
