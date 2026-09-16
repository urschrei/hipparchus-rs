#include <stdio.h>
#include <stdint.h>
#include "algebra.c"
static uint64_t s = 0x9E3779B97F4A7C15ULL;
static double rnd(void) { s ^= s << 13; s ^= s >> 7; s ^= s << 17; return (s >> 11) * (1.0 / 9007199254740992.0); }
int main(void) {
   FILE *f = fopen("stereo.csv", "w");
   for (int i = 0; i < 2000; i++) {
      struct ltln a = {(rnd() * 2 - 1) * 89.0 * DEG2RAD, (rnd() * 2 - 1) * PI};
      struct vct3 p0, pe, back; struct plpt w;
      LatLongToDcos3(&a, &p0);
      /* point within ~60 degrees of p0 */
      pe.di = p0.di + (rnd() * 2 - 1); pe.dj = p0.dj + (rnd() * 2 - 1); pe.dk = p0.dk + (rnd() * 2 - 1);
      NormalizeDcos3(&pe);
      MapStereo(&p0, &pe, &w);
      UnMapStereo(&p0, &w, &back);
      fprintf(f, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
         p0.di, p0.dj, p0.dk, pe.di, pe.dj, pe.dk, w.est, w.nrt, back.di, back.dj, back.dk, ArcDist(&p0, &pe));
   }
   fclose(f);
   return 0;
}
