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
   /* decode first vertex of every dense segment, as VIEW.C does */
   FILE *fi = fopen("COAST0.IDX", "rb"), *fd = fopen("COAST.DAT", "rb");
   f = fopen("coast0_first.csv", "w");
   double c[4]; int32_t off; int16_t n, id; short v[2];
   while (fread(c, 8, 4, fi) == 4 && fread(&off, 4, 1, fi) && fread(&n, 2, 1, fi) && fread(&id, 2, 1, fi)) {
      struct vct3 ctr = {c[0], c[1], c[2]}, out; struct plpt w;
      fseek(fd, off, SEEK_SET); fread(v, 2, 2, fd);
      w.est = c[3] / LC_SCALE * v[0]; w.nrt = c[3] / LC_SCALE * v[1];
      UnMapStereo(&ctr, &w, &out);
      fprintf(f, "%d,%.17g,%.17g,%.17g\n", off, out.di, out.dj, out.dk);
   }
   fclose(f);
   return 0;
}
