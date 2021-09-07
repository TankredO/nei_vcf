#ifdef _WIN32
#define LIBRARY_API extern "C" __declspec(dllexport)
#else
#define LIBRARY_API extern "C"
#endif

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<errno.h>

// 1MB (-1B for \0) line buffer
// TODO: maybe reduce buffer size, or make this adaptive
char line[1048576];

// maximum sample name length
const int MAX_SAMPLE_NAME_LEN = 256;

// source https://stackoverflow.com/a/8514474/5665958
char* mystrsep(char** stringp, const char* delim) {
  char* start = *stringp;
  char* p;

  p = (start != NULL) ? strpbrk(start, delim) : NULL;

  if (p == NULL) {
    *stringp = NULL;
  } else {
    *p = '\0';
    *stringp = p + 1;
  }

  return start;
}

LIBRARY_API void vcf_stats(char* vcf_file_path, int* shape, int n_unused_columns) {
    shape[0] = 0; shape[1] = 0;

    // open file
    FILE *fp;
    fp = fopen (vcf_file_path,"r");
    if (fp == NULL) {
        fprintf(stderr, "Problem opening file %s, errno = %d\n", vcf_file_path, errno);
        exit(1);
    }
    
    // iterate over lines in file
    int c = 0;
    char delimiters[] = " \t";
    char* running;
    char* token;
    while (fgets(line, sizeof line, fp)) {
        c++;
        size_t len = strlen(line);
        if (len && (line[len - 1] != '\n')) {
            fprintf(stderr, "Line %d exceeds buffer size!", c);
            exit(1);
        }

        // empty line
        if (line[0] == '\n') {
            continue;
        }

        // definitions (header)
        if (line[0] == '#' && line[1] == '#') {
            continue;
        }

        // table header
        int tc = 0;
        if (line[0] == '#') {
            running = line;
            while(token = mystrsep(&running, delimiters)) {
                if (tc > (n_unused_columns - 1)) {
                    // printf("%s\n", token);
                    shape[1] += 1;
                }
                tc += 1;
            }
        } else { 
            shape[0] += 1;
        }
    }

    // close file
    fclose(fp);
}

inline void str_to_variants(char *variant_str, double* variants) {
    char* running = variant_str;
    char* token;

    int a = atoi(mystrsep(&running, ","));
    int b = atoi(mystrsep(&running, ","));
    int c = atoi(mystrsep(&running, ","));
    int d = atoi(mystrsep(&running, ","));
    int sum = a + b + c + d;

    if (sum == 0) {
        variants[0] = 0.0;
        variants[1] = 0.0;
        variants[2] = 0.0;
        variants[3] = 0.0;
    } else {
        variants[0] = ((double) a) / sum;
        variants[1] = ((double) b) / sum;
        variants[2] = ((double) c) / sum;
        variants[3] = ((double) d) / sum;
    }
}

LIBRARY_API void get_variants(
    char* vcf_file_path,
    double* variants,
    char* sample_names,
    int* shape,
    int n_unused_columns,
    int entry_idx
) {
    // open file
    FILE *fp;
    fp = fopen (vcf_file_path,"r");
    if (fp == NULL) {
        fprintf(stderr, "Problem opening file %s, errno = %d\n", vcf_file_path, errno);
        exit(1);
    }

    // iterate over lines in file
    int c = 0;
    int rc = 0; // row counter
    char delimiters[] = " \t";
    char* running;
    char* running2;
    char* token;
    char* token2;
    while (fgets(line, sizeof line, fp)) {
        c++;
        size_t len = strlen(line);
        if (len && (line[len - 1] != '\n')) {
            fprintf(stderr, "Line %d exceeds buffer size!", c);
            exit(1);
        }

        // empty line
        if (line[0] == '\n') continue;
        // definitions (header)
        if (line[0] == '#' && line[1] == '#') continue;

        // table header
        int cc = 0; // column counter
        int sc = 0; // sample counter
        if (line[0] == '#') {
            running = line;
            while(token = mystrsep(&running, delimiters)) {
                if (cc > (n_unused_columns - 1)) {
                    strcpy(sample_names + sc * (MAX_SAMPLE_NAME_LEN+1), token);
                    sc += 1;
                }
                cc += 1;
            }
        // body
        } else {
            running = line;
            while(token = mystrsep(&running, delimiters)) {
                if (cc > (n_unused_columns - 1)) {
                    running2 = token;

                    int idx = sc * 4 + rc * shape[1] * 4;
                    
                    // entry "." means missing data
                    if (running2[0] == '.') {
                        *(variants + idx) = 0.0;
                        *(variants + idx + 1) = 0.0;
                        *(variants + idx + 2) = 0.0;
                        *(variants + idx + 3) = 0.0;
                        mystrsep(&running2, "\t");
                    } else {
                        // skip to entry_idx
                        for (int e_i = 0; e_i < entry_idx; e_i++) {
                            mystrsep(&running2, ":");
                        }

                        token2 = mystrsep(&running2, ":");

                        str_to_variants(token2, variants + idx);
                    }
                    sc += 1;
                }
                cc++;
            }
            rc++;
        }
    }

    fclose(fp);
}
