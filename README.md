# SVD image compressor

[![CI](https://github.com/buklinfur/SVD-img-compressor/actions/workflows/ci.yml/badge.svg)](https://github.com/buklinfur/SVD-img-compressor/actions/workflows/ci.yml)

A lossy image codec built on the truncated singular value decomposition, written
in plain C11 with **no linear-algebra libraries**: both SVD solvers (an exact
one-sided Jacobi method and a randomized top-k method) are implemented from
scratch on top of the C standard library and `libm`. Images are read and written
with the single-header [stb](third_party/stb) libraries, which are vendored, so
a fresh clone builds with one command.

![original and reconstructions with k = 5, 20, 50](docs/grid.jpg)

*Rank-k approximations of a 400×469 photo, per RGB channel. File sizes:
k = 5 → 12.9 KiB, k = 20 → 51.7 KiB, k = 50 → 129.1 KiB (raw pixels: 549.6 KiB).*

- Own SVD: Householder QR with column pivoting + one-sided Jacobi rotations
  (the structure of LAPACK's `xGEJSV`), and a randomized range finder with power
  iterations (Halko–Martinsson–Tropp) for the top-k case
- Agrees with `numpy.linalg.svd` to ~1e-13 (relative to σ₁), checked in CI
- Real compressed format (`.svdc`) with quantized factors, `compress`, `decompress`,
  `info` and `bench` commands
- Any image size, grayscale / gray+alpha / RGB / RGBA, PNG, JPEG, BMP, TGA input
- Clean under `-Wall -Wextra -Wpedantic -Wconversion -Werror`, AddressSanitizer,
  UBSan and valgrind; CI builds and tests it on Linux (gcc, clang), macOS (clang)
  and Windows (MSVC)
- Optional OpenMP parallelism over channels

## How it works

A grayscale image is an m×n matrix A. Its singular value decomposition

    A = U Σ Vᵀ = σ₁ u₁v₁ᵀ + σ₂ u₂v₂ᵀ + … ,   σ₁ ≥ σ₂ ≥ … ≥ 0

writes it as a sum of rank-one "layers" ordered by importance. Keeping only the
first k terms gives A_k, and by the **Eckart–Young theorem** no rank-k matrix is
closer to A: ‖A − A_k‖²_F = σ²_{k+1} + … + σ²_r. Storing A_k takes
k·(m + n + 1) numbers instead of m·n, which is the whole idea of the codec: each
channel is factorized separately, and only U_k, Σ_k and V_k are written to disk.

### SVD solver 1: one-sided Jacobi (exact, `--method jacobi`)

Hestenes' method orthogonalizes the columns of A by plane rotations; when all
columns are mutually orthogonal, their norms are the singular values. To make
it fast it is run on a preconditioned matrix (Drmač & Veselić 2008):

1. Householder QR with column pivoting, A·P = Q·R (|R_jj| non-increasing).
2. One-sided Jacobi on X = Rᵀ: for each column pair (p, q) the rotation angle
   that zeroes ⟨x_p, x_q⟩ comes from t² + 2ζt − 1 = 0, ζ = (‖x_q‖² − ‖x_p‖²) / (2⟨x_p, x_q⟩),
   taking the smaller root (Rutishauser), so the rotation never exceeds 45°.
3. **Stopping rule:** a pair is considered orthogonal when
   |⟨x_p, x_q⟩| ≤ m·ε·‖x_p‖·‖x_q‖ (ε = machine epsilon); iteration stops after
   the first sweep in which every pair passes (a safety cap of 40 sweeps;
   an 868×1541 photo channel needs 11).
4. U = Q·V_x and V = P·U_x.

The method is backward stable and never forms AᵀA, so it does not square the
condition number. Cost is O(mn²) for the QR plus O(n³) per sweep. Rotating Rᵀ
instead of A, column pivoting and a vectorizable dot product took the exact SVD
of an 868×1541 channel from 20.7 s to 7.8 s.

### SVD solver 2: randomized top-k (default, `--method randomized`)

For compression only the leading k triplets matter, so the default solver
(Halko, Martinsson & Tropp 2011, algorithms 4.4 and 5.1) works on a sketch:

1. Y = A·Ω with a Gaussian n×(k + p) matrix Ω (p = 10 oversampling columns;
   seeded splitmix64 + Box–Muller, so results are reproducible).
2. q = 2 power iterations Y ← A·(Aᵀ·Y), re-orthonormalizing after every product
   (Gram–Schmidt applied twice) so that small singular directions survive rounding.
3. B = Qᵀ·A is only (k + p)×n; its exact SVD is computed with the Jacobi solver
   above, and U = Q·U_B.

Cost is O(mnk). If k + p reaches min(m, n) the solver simply falls back to the
exact Jacobi SVD. Its error stays within a few percent of the Eckart–Young
optimum (0.2% on a real image channel, whose spectrum decays quickly).

### Choosing k

`-k N` keeps N triplets per channel. `--energy F` instead picks, per channel,
the smallest k with σ₁² + … + σ_k² ≥ F·‖A‖²_F; with the randomized solver the
sketch size doubles until the target is met. Note that for raw pixel intensities
the first triplet carries the mean brightness, so even k = 1 holds 70–80% of the
energy: useful targets are 0.99 and above (see [limitations](#limitations)).

## Building

Requirements: a C11 compiler and CMake ≥ 3.14 (or just GNU make).

```sh
git clone https://github.com/buklinfur/SVD-img-compressor.git
cd SVD-img-compressor
cmake -B build && cmake --build build      # produces build/svdc
ctest --test-dir build                    # unit, CLI and numpy cross-check tests
```

Without CMake: `make` (builds `build-make/svdc`), `make test`, `make OPENMP=1`.

| CMake option | default | effect |
|---|---|---|
| `SVDC_OPENMP` | `ON` | compress color channels in parallel if OpenMP is found |
| `SVDC_WERROR` | `OFF` | treat warnings as errors (used in CI) |
| `SVDC_SANITIZE` | `OFF` | AddressSanitizer + UndefinedBehaviorSanitizer |
| `SVDC_NATIVE` | `OFF` | `-march=native` (≈20% faster Jacobi on AVX2 machines) |
| `SVDC_BUILD_TESTS` | `ON` | build the tests |

## Usage

```text
$ svdc compress examples/hopper.png hopper.svdc -k 50
examples/hopper.png: 400x469, 3 channels -> hopper.svdc
  k = 50, 50, 50 (randomized SVD, 8-bit factors)
  size 129.1 KiB, 4.3x smaller than raw pixels (549.6 KiB), 2.32x vs input file (299.8 KiB)
  PSNR 27.57 dB, time 0.114 s

$ svdc decompress hopper.svdc hopper_k50.png
hopper.svdc -> hopper_k50.png (400x469, 3 channels), 0.033 s

$ svdc compress examples/hopper.png hopper.svdc --energy 0.99
  k = 46, 56, 46 (randomized SVD, 8-bit factors)
  ...

$ svdc bench examples/hopper.png --ks 5,20,50,100 --csv results.csv
examples/hopper.png: 400x469, 3 channels, raw 549.6 KiB, file 299.8 KiB, randomized SVD, 8-bit factors
     k         size  ratio(raw)  ratio(file)   PSNR dB  compress s  decompress s
     5     12.9 KiB       42.5x       23.17x     17.06       0.019         0.004
    20     51.7 KiB       10.6x        5.80x     22.44       0.035         0.008
    50    129.1 KiB        4.3x        2.32x     27.57       0.094         0.028
   100    258.1 KiB        2.1x        1.16x     32.86       0.205         0.046
```

The reported PSNR is always measured on what a decoder will produce, i.e. after
quantization. Other options: `--method jacobi|randomized`, `--bits 8|16|32`,
`--power-iters`, `--oversample`, `--seed`, `--threads`, `-q`; `svdc --help`
lists them all. `bench --save prefix` also writes the reconstructions as PNG.
Constant channels (such as a fully opaque alpha channel) are stored as a single
byte. Exit status is 0 on success, 1 on runtime errors and 2 on usage errors.

## Results

All numbers below are produced by `svdc bench` and
[`scripts/make_figures.py`](scripts/make_figures.py), which regenerates the
figures and tables in this README. Machine: 2 vCPU Intel Xeon @ 2.8 GHz,
gcc 13, `-O3`, no `-march=native`.

**Size and quality vs k** (`examples/hopper.png`, 400×469 RGB, 8-bit factors,
randomized SVD; raw pixels are 549.6 KiB):

| k | file size | vs raw | bits/pixel | PSNR, dB | compress, s |
|--:|--:|--:|--:|--:|--:|
| 1 | 2.6 KiB | 210.5x | 0.11 | 13.35 | 0.015 |
| 5 | 12.9 KiB | 42.5x | 0.56 | 17.07 | 0.019 |
| 10 | 25.8 KiB | 21.3x | 1.13 | 19.53 | 0.022 |
| 20 | 51.7 KiB | 10.6x | 2.26 | 22.44 | 0.034 |
| 50 | 129.1 KiB | 4.3x | 5.64 | 27.57 | 0.084 |
| 100 | 258.1 KiB | 2.1x | 11.27 | 32.86 | 0.205 |
| 200 | 516.2 KiB | 1.1x | 22.54 | 41.32 | 0.740 |

![PSNR vs bits per pixel for svdc and JPEG](docs/rate_distortion.png)

**Honest comparison with JPEG.** At equal size JPEG is far better: it reaches
~29 dB at 0.5 bits per pixel, while the SVD codec needs ~6 bits per pixel for
the same quality. This is expected and is the main lesson of the project: a
truncated SVD is the *optimal low-rank* approximation, but natural images are
not low-rank. Their singular values decay slowly (below), and fine detail costs
whole rank-one layers spanning the full image, while JPEG spends bits locally
(8×8 DCT blocks), exploits human perception (chroma subsampling, quantization
tables) and entropy-codes the result. Here 8-bit factors beat 16-bit ones at
every size, because quantization noise is much smaller than the truncation error.

![singular values and captured energy per channel](docs/spectrum.png)

**Speed** (868×1541 RGB photo, `examples/img.jpg`, k = 50 unless noted):

| solver | 1 thread | 2 threads (OpenMP over channels) |
|---|--:|--:|
| randomized, k = 50 | 0.94 s | 0.45 s |
| randomized, k = 200 | 4.16 s | 2.63 s |
| Jacobi (exact, all 868 triplets) | 21.0 s | 14.5 s |

Three channels on two cores bound the channel-level speedup at 1.5× for the
compute-bound Jacobi solver. The randomized solver loses at most 0.15 dB of
PSNR against the exact one (hopper.png, k = 5…200; the exact SVD of that image
takes 1.0 s, the randomized one 0.02–0.65 s depending on k).

**Accuracy** (`tests/numpy_check.py`, run in CI): the Jacobi SVD matches
`numpy.linalg.svd` with max |σᵢ − σᵢ^np| / σ₁ ≤ 6e-14,
‖A − UΣVᵀ‖/‖A‖ ≤ 1.2e-13 and ‖UᵀU − I‖, ‖VᵀV − I‖ ≤ 1.7e-13, including a matrix
with condition number 1e12 and a rank-deficient one. The randomized top-20
error is at most 1.4% above the Eckart–Young optimum (Gaussian noise, the worst
case for sketching) and 0.2% above it on an image channel.

## File format (`.svdc`, version 1)

All integers are little-endian, floats are IEEE-754 binary32.

| offset | size | field |
|--:|--:|---|
| 0 | 4 | magic `SVDC` |
| 4 | 1 | version (1) |
| 5 | 1 | bits per factor value: 8, 16 or 32 |
| 6 | 1 | number of channels, 1–4 |
| 7 | 1 | reserved (0) |
| 8 | 4 | width (u32) |
| 12 | 4 | height (u32) |
| 16 | … | one block per channel |

Channel block:

| field | size | notes |
|---|---|---|
| type | u8 | 0 = constant channel, 1 = low rank |
| value | u8 | constant channels only |
| k | u32 | low rank: number of triplets, 1 ≤ k ≤ min(width, height) |
| σ₁ … σ_k | k × f32 | singular values |
| U columns | k × (scale + height values) | column-by-column |
| V columns | k × (scale + width values) | column-by-column |

With 8 or 16 bits, each column stores an f32 scale (its max |entry|) followed by
signed integers q, decoded as q / 127 · scale or q / 32767 · scale; with 32 bits
the values are raw f32 and there is no scale. A low-rank channel therefore costs
about k·(width + height) bytes at 8 bits. The parser validates every field and
the remaining length before allocating, and rejects truncated files and
trailing data (the tests try every prefix and every single-byte corruption of a
file).

## Project layout

```text
src/matrix.[ch]    dense row-major matrix, products, Gram-Schmidt
src/svd.[ch]       Jacobi and randomized SVD, energy-based rank selection
src/image.[ch]     image I/O (stb), channel <-> matrix conversion
src/codec.[ch]     per-channel compression, .svdc encoder and decoder
src/metrics.[ch]   MSE / PSNR, wall-clock timer
src/main.c         command line interface
src/stb_impl.c     the only translation unit that compiles stb
tests/             unit tests, CLI tests, numpy cross-check
scripts/           README figure generator
third_party/stb/   vendored stb_image / stb_image_write (MIT / public domain)
```

## Limitations

- **It is not competitive with JPEG/WebP** (see above); it is a study of
  low-rank approximation, not a production codec.
- `--energy` is measured on raw intensities, where the mean dominates; targets
  below ~0.99 give very small k.
- No entropy coding and no color transform. Converting to YCbCr with a smaller
  k for chroma, subtracting the per-column mean, and compressing the quantized
  factors with an entropy coder would each shrink files noticeably.
- The whole image is held in memory as doubles (8 bytes per pixel per channel);
  the exact Jacobi solver also needs O(min(m, n)²) memory and is slow for
  images larger than a few megapixels.
- 8-bit factors cap quality at about 50 dB; use `--bits 16` if that matters.
- Only 8-bit images: 16-bit PNGs are converted to 8 bits on load.

## History

This started as my first-year C course project (2024), which called LAPACKE
and only worked for square PNG images. That version is kept unchanged on the
[`legacy`](https://github.com/buklinfur/SVD-img-compressor/tree/legacy) branch.

## License

MIT, see [LICENSE](LICENSE). stb is MIT / public domain
([third_party/stb/LICENSE](third_party/stb/LICENSE)). Image credits are in
[examples/README.md](examples/README.md).
