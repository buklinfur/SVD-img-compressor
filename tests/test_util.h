#ifndef SVDC_TEST_UTIL_H
#define SVDC_TEST_UTIL_H

// minimal assertion helpers for the ctest executables

#include <math.h>
#include <stdio.h>

static int testFailures = 0;

#define CHECK(cond)                                                           \
    do {                                                                      \
        if (!(cond)) {                                                        \
            fprintf (stderr, "%s:%d: CHECK failed: %s\n", __FILE__, __LINE__, \
                     #cond);                                                  \
            testFailures++;                                                   \
        }                                                                     \
    } while (0)

#define CHECK_NEAR(actual, expected, tol)                                       \
    do {                                                                        \
        double a_ = (double) (actual);                                          \
        double e_ = (double) (expected);                                        \
        if (!(fabs (a_ - e_) <= (tol))) {                                       \
            fprintf (stderr, "%s:%d: CHECK_NEAR failed: %s = %.17g, want %.17g" \
                     " (tol %.3g)\n", __FILE__, __LINE__, #actual, a_, e_,     \
                     (double) (tol));                                           \
            testFailures++;                                                     \
        }                                                                       \
    } while (0)

#define CHECK_LESS(actual, bound)                                               \
    do {                                                                        \
        double a_ = (double) (actual);                                          \
        double b_ = (double) (bound);                                           \
        if (!(a_ < b_)) {                                                       \
            fprintf (stderr, "%s:%d: CHECK_LESS failed: %s = %.6g, want < %.6g\n", \
                     __FILE__, __LINE__, #actual, a_, b_);                     \
            testFailures++;                                                     \
        }                                                                       \
    } while (0)

#define RUN_TEST(fn)                         \
    do {                                     \
        int before_ = testFailures;          \
        fn ();                               \
        printf ("%-40s %s\n", #fn,           \
                testFailures == before_ ? "ok" : "FAILED"); \
    } while (0)

#define TEST_RESULT() (testFailures == 0 ? 0 : 1)

#endif
