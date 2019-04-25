#ifndef MAIN_H
#define MAIN_H

#define	CUBE(x) (x*x*x)
#define MAX_STRING_LEN 1024

#include <stdint.h>

// Structs //

// Proto-Types //


// Macros //
// If the value of `status` is less than 0 (i.e., a failure state), `status` is returned.
// Also prints the passed error message to stdout.
#define CHECK_STATUS_AND_RETURN_ON_FAIL(status, return_value, ...) \
    do {                                                           \
        if(status < 0) {                                           \
            fprintf(stderr, __VA_ARGS__);                          \
            return status;                                         \
        }                                                          \
  } while (0)

// As above, but no printing.
#define CHECK_STATUS_AND_RETURN_ON_FAIL_NO_PRINT(status, return_value) \
    do {                                                               \
        if(status < 0) {                                               \
            return status;                                             \
        }                                                              \
  } while (0)
#endif
