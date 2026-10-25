#!/usr/bin/env python3
"""Cross-check the C SVD implementations against numpy.linalg.svd.

usage: numpy_check.py --tool build/svd_tool [--image examples/test2.png]

Exits with status 1 if any check exceeds its tolerance.
"""

import argparse
import os
import subprocess
import sys
import tempfile

import numpy as np


def run_tool(tool, method, k, a, workdir):
    src = os.path.join(workdir, "in.bin")
    dst = os.path.join(workdir, "out.bin")
    with open(src, "wb") as f:
        np.array(a.shape, dtype=np.int32).tofile(f)
        np.ascontiguousarray(a, dtype=np.float64).tofile(f)
    subprocess.run([tool, method, str(k), src, dst], check=True)
    with open(dst, "rb") as f:
        rank, rows, cols = np.fromfile(f, dtype=np.int32, count=3)
        s = np.fromfile(f, dtype=np.float64, count=rank)
        u = np.fromfile(f, dtype=np.float64, count=rows * rank).reshape(rows, rank)
        v = np.fromfile(f, dtype=np.float64, count=cols * rank).reshape(cols, rank)
    return s, u, v


def matrix_with_spectrum(rng, m, n, s):
    q1, _ = np.linalg.qr(rng.standard_normal((m, len(s))))
    q2, _ = np.linalg.qr(rng.standard_normal((n, len(s))))
    return (q1 * s) @ q2.T


def cases(rng, image):
    yield "gaussian 200x120", rng.standard_normal((200, 120))
    yield "gaussian 90x300", rng.standard_normal((90, 300))
    yield "graded 150x100 (cond 1e12)", matrix_with_spectrum(
        rng, 150, 100, np.logspace(0, -12, 100))
    yield "rank 7 of 80x60", matrix_with_spectrum(rng, 80, 60, np.arange(7, 0, -1.0))
    yield "geometric decay 400x300", matrix_with_spectrum(
        rng, 400, 300, 100 * 0.97 ** np.arange(300))
    if image:
        from PIL import Image  # only needed for the optional image case
        pixels = np.asarray(Image.open(image).convert("RGB"), dtype=np.float64)
        yield "image %s (red channel %dx%d)" % (
            os.path.basename(image), pixels.shape[0], pixels.shape[1]), pixels[..., 0]


def check(name, value, tol, failures):
    ok = value <= tol
    print("    %-34s %.3e  (tol %.0e) %s" % (name, value, tol, "ok" if ok else "FAIL"))
    if not ok:
        failures.append(name)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--tool", required=True)
    parser.add_argument("--image")
    parser.add_argument("--k", type=int, default=20)
    args = parser.parse_args()

    rng = np.random.default_rng(2024)
    failures = []
    with tempfile.TemporaryDirectory() as workdir:
        for title, a in cases(rng, args.image):
            print(title)
            s_np = np.linalg.svd(a, compute_uv=False)
            norm = np.linalg.norm(a)

            # exact Jacobi SVD: compare everything
            s, u, v = run_tool(args.tool, "jacobi", 0, a, workdir)
            r = len(s)
            nonzero = s > 0
            check("jacobi |s - s_np|_max / s_1", np.abs(s - s_np[:r]).max() / s_np[0], 1e-12, failures)
            check("jacobi ||A - U S V^T|| / ||A||",
                  np.linalg.norm(a - (u * s) @ v.T) / norm, 1e-12, failures)
            un, vn = u[:, nonzero], v[:, nonzero]
            check("jacobi ||U^T U - I||_max",
                  np.abs(un.T @ un - np.eye(un.shape[1])).max(), 1e-11, failures)
            check("jacobi ||V^T V - I||_max",
                  np.abs(vn.T @ vn - np.eye(vn.shape[1])).max(), 1e-11, failures)

            # randomized top-k: compare with the optimal rank-k error (Eckart-Young)
            k = min(args.k, min(a.shape))
            s, u, v = run_tool(args.tool, "randomized", k, a, workdir)
            optimal = np.sqrt(np.sum(s_np[k:] ** 2))
            err = np.linalg.norm(a - (u * s) @ v.T)
            ratio = err / optimal if optimal > 1e-12 * norm else err / norm
            check("randomized k=%d error / optimal - 1" % k, max(ratio - 1.0, 0.0), 0.05, failures)
            # individual singular values are only approximate for flat spectra
            # (Gaussian noise is the worst case); what matters is the error above
            check("randomized k=%d |s - s_np|_max / s_1" % k,
                  np.abs(s - s_np[:k]).max() / s_np[0], 1e-1, failures)

    if failures:
        print("\n%d check(s) failed" % len(failures))
        return 1
    print("\nall checks passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
