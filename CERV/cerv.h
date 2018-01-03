#ifndef CERV_H
#define CERV_H

#ifdef __cplusplus
extern "C" {
#endif
void cerv_encode(int** img, int NR, int NC, const char* filename);
void cerv_decode(int** img, int NR, int NC, const char* filename);
#ifdef __cplusplus
}
#endif

#endif
