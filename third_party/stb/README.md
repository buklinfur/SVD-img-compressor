# stb (vendored)

Single-header image I/O libraries by Sean Barrett and contributors.

| File                | Version | Purpose                         |
|---------------------|---------|---------------------------------|
| `stb_image.h`       | v2.30   | decoding PNG, JPEG, BMP, TGA... |
| `stb_image_write.h` | v1.16   | encoding PNG, JPEG, BMP, TGA    |

- Upstream: https://github.com/nothings/stb
- Snapshot: commit `003501f6d7a8668ce8b6d875e34dd97ad0d5d5b5`
- License: dual MIT / public domain (Unlicense), see [`LICENSE`](LICENSE)
  and the notice at the end of each header.

The headers are unmodified. Their implementation is compiled once, in
`src/stb_impl.c`, with compiler warnings disabled for that translation unit
only, so the project's own strict warning flags do not apply to third-party code.
