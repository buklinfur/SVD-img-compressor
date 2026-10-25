#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "codec.h"
#include "image.h"
#include "metrics.h"
#include "svd.h"

#ifndef SVDC_VERSION_STRING
#define SVDC_VERSION_STRING "1.0.0"
#endif

#define EXIT_USAGE 2
#define MAX_KS 64

static void printUsage (FILE* f) {
    fprintf (f,
        "svdc " SVDC_VERSION_STRING " - image compression with a truncated SVD\n"
        "\n"
        "usage:\n"
        "  svdc compress   <image> <output.svdc> [options]\n"
        "  svdc decompress <input.svdc> <image.png|jpg|bmp|tga> [-q]\n"
        "  svdc info       <input.svdc>\n"
        "  svdc bench      <image> [--ks 5,20,50] [--csv file.csv] [--save prefix] [options]\n"
        "\n"
        "compression options:\n"
        "  -k N               singular triplets kept per channel (default 50)\n"
        "  --energy F         instead of -k: keep the smallest k with F of the energy\n"
        "                     (sum of squared singular values), 0 < F <= 1\n"
        "  --method M         randomized (default) or jacobi (exact, slower)\n"
        "  --bits B           precision of stored factors: 8 (default), 16 or 32\n"
        "  --power-iters N    randomized: power iterations (default 2)\n"
        "  --oversample N     randomized: extra sketch columns (default 10)\n"
        "  --seed N           randomized: seed of the random sketch (default 1)\n"
        "  --threads N        threads for the per-channel loop (OpenMP builds)\n"
        "  -q, --quiet        print nothing on success\n");
}

// ---------------------------------------------------------------------------
// argument parsing

typedef struct Args {
    const char* positional[2];
    int positionalCount;
    CompressOptions opt;
    int quiet;
    int kGiven;
    int ks[MAX_KS];
    int ksCount;
    const char* csvPath;
    const char* savePrefix;
} Args;

static int parseInt (const char* text, long lo, long hi, int* out) {
    char* end = NULL;
    errno = 0;
    long v = strtol (text, &end, 10);
    if (errno != 0 || end == text || *end != '\0' || v < lo || v > hi)
        return -1;
    *out = (int) v;
    return 0;
}

static int parseDouble (const char* text, double* out) {
    char* end = NULL;
    errno = 0;
    double v = strtod (text, &end);
    if (errno != 0 || end == text || *end != '\0' || !isfinite (v))
        return -1;
    *out = v;
    return 0;
}

static int parseKs (const char* text, Args* args) {
    args->ksCount = 0;
    const char* p = text;
    while (*p != '\0') {
        char item[32];
        size_t len = strcspn (p, ",");
        if (len == 0 || len >= sizeof item || args->ksCount >= MAX_KS)
            return -1;
        memcpy (item, p, len);
        item[len] = '\0';
        if (parseInt (item, 1, 1L << 24, &args->ks[args->ksCount]) != 0)
            return -1;
        args->ksCount++;
        p += len;
        if (*p == ',')
            p++;
    }
    return args->ksCount > 0 ? 0 : -1;
}

static int parseArgs (int argc, char** argv, Args* args) {
    memset (args, 0, sizeof *args);
    args->opt = defaultCompressOptions ();

    for (int i = 2; i < argc; i++) {
        const char* a = argv[i];
        const char* value = i + 1 < argc ? argv[i + 1] : NULL;
        int ok = 0;

        if (strcmp (a, "-q") == 0 || strcmp (a, "--quiet") == 0) {
            args->quiet = 1;
            continue;
        }
        if (a[0] != '-' || a[1] == '\0') {
            if (args->positionalCount >= 2) {
                fprintf (stderr, "svdc: unexpected argument '%s'\n", a);
                return -1;
            }
            args->positional[args->positionalCount++] = a;
            continue;
        }
        if (value == NULL) {
            fprintf (stderr, "svdc: option %s needs a value\n", a);
            return -1;
        }
        i++;

        if (strcmp (a, "-k") == 0) {
            ok = parseInt (value, 1, 1L << 24, &args->opt.k) == 0;
            args->kGiven = 1;
        } else if (strcmp (a, "--energy") == 0) {
            ok = parseDouble (value, &args->opt.energy) == 0 && args->opt.energy > 0.0 &&
                 args->opt.energy <= 1.0;
        } else if (strcmp (a, "--method") == 0) {
            ok = 1;
            if (strcmp (value, "jacobi") == 0)
                args->opt.svd.method = SVD_JACOBI;
            else if (strcmp (value, "randomized") == 0)
                args->opt.svd.method = SVD_RANDOMIZED;
            else
                ok = 0;
        } else if (strcmp (a, "--bits") == 0) {
            ok = parseInt (value, 8, 32, &args->opt.bits) == 0 &&
                 (args->opt.bits == 8 || args->opt.bits == 16 || args->opt.bits == 32);
        } else if (strcmp (a, "--power-iters") == 0) {
            ok = parseInt (value, 0, 20, &args->opt.svd.powerIters) == 0;
        } else if (strcmp (a, "--oversample") == 0) {
            ok = parseInt (value, 0, 1000, &args->opt.svd.oversampling) == 0;
        } else if (strcmp (a, "--seed") == 0) {
            int seed = 0;
            ok = parseInt (value, 0, 2147483647L, &seed) == 0;
            args->opt.svd.seed = (uint64_t) seed;
        } else if (strcmp (a, "--threads") == 0) {
            ok = parseInt (value, 1, 1024, &args->opt.threads) == 0;
        } else if (strcmp (a, "--ks") == 0) {
            ok = parseKs (value, args) == 0;
        } else if (strcmp (a, "--csv") == 0) {
            args->csvPath = value;
            ok = 1;
        } else if (strcmp (a, "--save") == 0) {
            args->savePrefix = value;
            ok = 1;
        } else {
            fprintf (stderr, "svdc: unknown option '%s'\n", a);
            return -1;
        }
        if (!ok) {
            fprintf (stderr, "svdc: invalid value '%s' for %s\n", value, a);
            return -1;
        }
    }
    if (args->kGiven && args->opt.energy > 0.0) {
        fprintf (stderr, "svdc: use either -k or --energy, not both\n");
        return -1;
    }
    return 0;
}

// ---------------------------------------------------------------------------
// helpers

static long fileSize (const char* filename) {
    FILE* f = fopen (filename, "rb");
    if (f == NULL)
        return -1;
    long size = -1;
    if (fseek (f, 0, SEEK_END) == 0)
        size = ftell (f);
    fclose (f);
    return size;
}

static void formatBytes (double bytes, char* out, size_t size) {
    if (bytes < 1024.0)
        snprintf (out, size, "%.0f B", bytes);
    else if (bytes < 1024.0 * 1024.0)
        snprintf (out, size, "%.1f KiB", bytes / 1024.0);
    else
        snprintf (out, size, "%.2f MiB", bytes / (1024.0 * 1024.0));
}

static const char* methodName (SvdMethod m) {
    return m == SVD_JACOBI ? "jacobi" : "randomized";
}

static size_t rawSize (const Image* img) {
    return (size_t) img->width * (size_t) img->height * (size_t) img->channels;
}

// "k = 20, 20, 20" or "k = 20, const" for the summary line
static void describeRanks (const Compressed* c, char* out, size_t size) {
    size_t pos = (size_t) snprintf (out, size, "k =");
    for (int ch = 0; ch < c->channels && pos < size; ch++) {
        const CompressedChannel* cc = &c->channel[ch];
        if (cc->constant || cc->svd == NULL)
            pos += (size_t) snprintf (out + pos, size - pos, "%s const", ch ? "," : "");
        else
            pos += (size_t) snprintf (out + pos, size - pos, "%s %d", ch ? "," : "", cc->svd->rank);
    }
}

static int writeBuffer (const char* filename, const unsigned char* buf, size_t size) {
    FILE* f = fopen (filename, "wb");
    if (f == NULL)
        return -1;
    size_t done = fwrite (buf, 1, size, f);
    int closed = fclose (f);
    return done == size && closed == 0 ? 0 : -1;
}

static void printPsnr (FILE* f, double psnr) {
    if (isinf (psnr))
        fprintf (f, "inf");
    else
        fprintf (f, "%.2f", psnr);
}

// ---------------------------------------------------------------------------
// commands

static int cmdCompress (const Args* args) {
    if (args->positionalCount != 2) {
        fprintf (stderr, "usage: svdc compress <image> <output.svdc> [options]\n");
        return EXIT_USAGE;
    }
    const char* input = args->positional[0];
    const char* output = args->positional[1];

    Image* img = readImage (input);
    if (img == NULL) {
        fprintf (stderr, "svdc: cannot read '%s': %s\n", input, imageFailureReason ());
        return EXIT_FAILURE;
    }

    double t0 = wallTime ();
    Compressed* c = compressImage (img, &args->opt);
    unsigned char* buf = NULL;
    size_t size = 0;
    int status = c != NULL ? encodeCompressed (c, &buf, &size) : SVDC_ERR_ALLOC;
    double elapsed = wallTime () - t0;
    if (status == SVDC_OK && writeBuffer (output, buf, size) != 0)
        status = SVDC_ERR_IO;
    if (status != SVDC_OK) {
        fprintf (stderr, "svdc: cannot compress '%s' to '%s': %s\n", input, output,
                 codecErrorString (status));
        free (buf);
        freeCompressed (c);
        freeImage (img);
        return EXIT_FAILURE;
    }

    if (!args->quiet) {
        // measure quality on exactly what a decoder will see (quantized factors)
        Compressed* decoded = NULL;
        Image* restored = NULL;
        double psnr = (double) NAN;
        if (decodeCompressed (buf, size, &decoded) == SVDC_OK)
            restored = decompressImage (decoded);
        if (restored != NULL)
            psnr = imagePsnr (img, restored);

        char ranks[96], sizeText[32], rawText[32], fileText[32];
        describeRanks (c, ranks, sizeof ranks);
        formatBytes ((double) size, sizeText, sizeof sizeText);
        formatBytes ((double) rawSize (img), rawText, sizeof rawText);
        long inputSize = fileSize (input);
        formatBytes ((double) inputSize, fileText, sizeof fileText);

        printf ("%s: %dx%d, %d channel%s -> %s\n", input, img->width, img->height,
                img->channels, img->channels == 1 ? "" : "s", output);
        printf ("  %s (%s SVD, %d-bit factors)\n", ranks, methodName (args->opt.svd.method), c->bits);
        printf ("  size %s, %.1fx smaller than raw pixels (%s), %.2fx vs input file (%s)\n",
                sizeText, (double) rawSize (img) / (double) size, rawText,
                inputSize > 0 ? (double) inputSize / (double) size : 0.0, fileText);
        printf ("  PSNR ");
        printPsnr (stdout, psnr);
        printf (" dB, time %.3f s\n", elapsed);

        freeImage (restored);
        freeCompressed (decoded);
    }

    free (buf);
    freeCompressed (c);
    freeImage (img);
    return EXIT_SUCCESS;
}

static int cmdDecompress (const Args* args) {
    if (args->positionalCount != 2) {
        fprintf (stderr, "usage: svdc decompress <input.svdc> <image.png|jpg|bmp|tga>\n");
        return EXIT_USAGE;
    }
    const char* input = args->positional[0];
    const char* output = args->positional[1];
    if (!isWritableImageName (output)) {
        fprintf (stderr, "svdc: '%s': output must end in .png, .jpg, .jpeg, .bmp or .tga\n", output);
        return EXIT_USAGE;
    }

    double t0 = wallTime ();
    Compressed* c = NULL;
    int status = readCompressedFile (input, &c, NULL);
    if (status != SVDC_OK) {
        fprintf (stderr, "svdc: cannot read '%s': %s\n", input, codecErrorString (status));
        return EXIT_FAILURE;
    }
    Image* img = decompressImage (c);
    double elapsed = wallTime () - t0;
    if (img == NULL || writeImage (img, output) != 0) {
        fprintf (stderr, "svdc: cannot write '%s'\n", output);
        freeImage (img);
        freeCompressed (c);
        return EXIT_FAILURE;
    }
    if (!args->quiet)
        printf ("%s -> %s (%dx%d, %d channel%s), %.3f s\n", input, output, img->width,
                img->height, img->channels, img->channels == 1 ? "" : "s", elapsed);
    freeImage (img);
    freeCompressed (c);
    return EXIT_SUCCESS;
}

static int cmdInfo (const Args* args) {
    if (args->positionalCount != 1) {
        fprintf (stderr, "usage: svdc info <input.svdc>\n");
        return EXIT_USAGE;
    }
    Compressed* c = NULL;
    size_t size = 0;
    int status = readCompressedFile (args->positional[0], &c, &size);
    if (status != SVDC_OK) {
        fprintf (stderr, "svdc: cannot read '%s': %s\n", args->positional[0], codecErrorString (status));
        return EXIT_FAILURE;
    }
    char sizeText[32], ranks[96];
    formatBytes ((double) size, sizeText, sizeof sizeText);
    describeRanks (c, ranks, sizeof ranks);
    size_t raw = (size_t) c->width * (size_t) c->height * (size_t) c->channels;
    printf ("%s: svdc v1, %dx%d, %d channel%s, %d-bit factors\n", args->positional[0],
            c->width, c->height, c->channels, c->channels == 1 ? "" : "s", c->bits);
    printf ("  %s\n  size %s, %.1fx smaller than raw pixels\n", ranks, sizeText,
            (double) raw / (double) size);
    freeCompressed (c);
    return EXIT_SUCCESS;
}

// deep copy of a compressed image keeping at most k triplets per channel
static Compressed* truncatedCopy (const Compressed* src, int k) {
    Compressed* c = (Compressed*) calloc (1, sizeof (Compressed));
    if (c == NULL)
        return NULL;
    *c = *src;
    for (int ch = 0; ch < 4; ch++)
        c->channel[ch].svd = NULL;
    for (int ch = 0; ch < src->channels; ch++) {
        const Svd* s = src->channel[ch].svd;
        if (s == NULL)
            continue;
        int kk = k < s->rank ? k : s->rank;
        Svd* d = makeSvd (s->rows, s->cols, kk);
        if (d == NULL) {
            freeCompressed (c);
            return NULL;
        }
        memcpy (d->s, s->s, (size_t) kk * sizeof (double));
        for (int i = 0; i < s->rows; i++)
            memcpy (d->u->data + (size_t) i * (size_t) kk, s->u->data + (size_t) i * (size_t) s->rank,
                    (size_t) kk * sizeof (double));
        for (int i = 0; i < s->cols; i++)
            memcpy (d->v->data + (size_t) i * (size_t) kk, s->v->data + (size_t) i * (size_t) s->rank,
                    (size_t) kk * sizeof (double));
        c->channel[ch].svd = d;
    }
    return c;
}

static int cmdBench (Args* args) {
    if (args->positionalCount != 1) {
        fprintf (stderr, "usage: svdc bench <image> [--ks 5,20,50] [--csv file] [--save prefix] [options]\n");
        return EXIT_USAGE;
    }
    const char* input = args->positional[0];
    Image* img = readImage (input);
    if (img == NULL) {
        fprintf (stderr, "svdc: cannot read '%s': %s\n", input, imageFailureReason ());
        return EXIT_FAILURE;
    }
    if (args->ksCount == 0) {
        const int defaults[] = { 1, 5, 10, 20, 50, 100, 200 };
        for (size_t i = 0; i < sizeof defaults / sizeof defaults[0]; i++)
            args->ks[args->ksCount++] = defaults[i];
    }
    const int minDim = img->width < img->height ? img->width : img->height;
    const long inputSize = fileSize (input);
    const size_t raw = rawSize (img);

    FILE* csv = NULL;
    if (args->csvPath != NULL) {
        csv = fopen (args->csvPath, "w");
        if (csv == NULL) {
            fprintf (stderr, "svdc: cannot write '%s'\n", args->csvPath);
            freeImage (img);
            return EXIT_FAILURE;
        }
        fprintf (csv, "k,bytes,ratio_raw,ratio_file,psnr_db,compress_s,decompress_s\n");
    }

    // the exact Jacobi SVD does not depend on k: factorize once, truncate per row
    Compressed* full = NULL;
    double fullTime = 0.0;
    if (args->opt.svd.method == SVD_JACOBI) {
        CompressOptions o = args->opt;
        o.k = minDim;
        o.energy = 0.0;
        double t0 = wallTime ();
        full = compressImage (img, &o);
        fullTime = wallTime () - t0;
        if (full == NULL) {
            fprintf (stderr, "svdc: out of memory\n");
            if (csv != NULL)
                fclose (csv);
            freeImage (img);
            return EXIT_FAILURE;
        }
    }

    char rawText[32], fileText[32];
    formatBytes ((double) raw, rawText, sizeof rawText);
    formatBytes ((double) inputSize, fileText, sizeof fileText);
    printf ("%s: %dx%d, %d channel%s, raw %s, file %s, %s SVD, %d-bit factors\n", input,
            img->width, img->height, img->channels, img->channels == 1 ? "" : "s", rawText,
            fileText, methodName (args->opt.svd.method), args->opt.bits);
    printf ("%6s %12s %11s %12s %9s %11s %13s\n", "k", "size", "ratio(raw)", "ratio(file)",
            "PSNR dB", "compress s", "decompress s");

    int status = EXIT_SUCCESS;
    for (int i = 0; i < args->ksCount && status == EXIT_SUCCESS; i++) {
        const int k = args->ks[i];
        if (k > minDim)
            continue;

        double t0 = wallTime ();
        Compressed* c;
        if (full != NULL) {
            c = truncatedCopy (full, k);
        } else {
            CompressOptions o = args->opt;
            o.k = k;
            o.energy = 0.0;
            c = compressImage (img, &o);
        }
        unsigned char* buf = NULL;
        size_t size = 0;
        int enc = c != NULL ? encodeCompressed (c, &buf, &size) : SVDC_ERR_ALLOC;
        double tc = wallTime () - t0 + (full != NULL ? fullTime : 0.0);

        t0 = wallTime ();
        Compressed* decoded = NULL;
        Image* restored = NULL;
        if (enc == SVDC_OK && decodeCompressed (buf, size, &decoded) == SVDC_OK)
            restored = decompressImage (decoded);
        double td = wallTime () - t0;

        if (restored == NULL) {
            fprintf (stderr, "svdc: benchmark failed at k = %d\n", k);
            status = EXIT_FAILURE;
        } else {
            double psnr = imagePsnr (img, restored);
            char sizeText[32];
            formatBytes ((double) size, sizeText, sizeof sizeText);
            double ratioRaw = (double) raw / (double) size;
            double ratioFile = inputSize > 0 ? (double) inputSize / (double) size : 0.0;
            printf ("%6d %12s %10.1fx %11.2fx %9.2f %11.3f %13.3f\n", k, sizeText, ratioRaw,
                    ratioFile, psnr, tc, td);
            if (csv != NULL)
                fprintf (csv, "%d,%zu,%.4f,%.4f,%.4f,%.6f,%.6f\n", k, size, ratioRaw, ratioFile,
                         psnr, tc, td);
            if (args->savePrefix != NULL) {
                char path[1024];
                snprintf (path, sizeof path, "%s_k%d.png", args->savePrefix, k);
                if (writeImage (restored, path) != 0) {
                    fprintf (stderr, "svdc: cannot write '%s'\n", path);
                    status = EXIT_FAILURE;
                }
            }
        }
        freeImage (restored);
        freeCompressed (decoded);
        free (buf);
        freeCompressed (c);
    }
    if (full != NULL)
        printf ("(jacobi: one full SVD of %.3f s is shared by all rows and included in compress s)\n", fullTime);

    if (csv != NULL && fclose (csv) != 0) {
        fprintf (stderr, "svdc: cannot write '%s'\n", args->csvPath);
        status = EXIT_FAILURE;
    }
    freeCompressed (full);
    freeImage (img);
    return status;
}

int main (int argc, char** argv) {
    if (argc < 2) {
        printUsage (stderr);
        return EXIT_USAGE;
    }
    const char* command = argv[1];
    if (strcmp (command, "-h") == 0 || strcmp (command, "--help") == 0 || strcmp (command, "help") == 0) {
        printUsage (stdout);
        return EXIT_SUCCESS;
    }
    if (strcmp (command, "--version") == 0) {
        printf ("svdc " SVDC_VERSION_STRING "\n");
        return EXIT_SUCCESS;
    }

    Args args;
    if (parseArgs (argc, argv, &args) != 0)
        return EXIT_USAGE;

    if (strcmp (command, "compress") == 0)
        return cmdCompress (&args);
    if (strcmp (command, "decompress") == 0)
        return cmdDecompress (&args);
    if (strcmp (command, "info") == 0)
        return cmdInfo (&args);
    if (strcmp (command, "bench") == 0)
        return cmdBench (&args);

    fprintf (stderr, "svdc: unknown command '%s' (try 'svdc --help')\n", command);
    return EXIT_USAGE;
}
