#ifndef UTILS_H
#define UTILS_H
#endif

#ifndef MAXLENGTH
#define MAXLENGTH 1024
#endif

#ifdef NDEBUG
#define XASSERT(EXP, ...)                                do{} while(0)
#else
#define XASSERT(EXP, ...)                                              \
    do { if (!(EXP)) {                                                  \
            printf("Error in file: %s\tfunc: %s\tline: %d with expression `"#EXP"'\n", __FILE__, __FUNCTION__, __LINE__); \
            printf(__VA_ARGS__);                                        \
            printf("Bug in code: email Anne Hutter <ahutter@swin.edu.au>\n"); \
            fflush(stdout);                                             \
            exit(EXIT_FAILURE);                                         \
        } \
    } while (0)
#endif

int file_exist(char *filename);
