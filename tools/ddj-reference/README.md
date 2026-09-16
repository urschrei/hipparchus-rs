# DDJ 1992 reference harness

`ref.c` generates `tests/data/stereo.csv` by calling the original `ALGEBRA.C`
from Lukatela & Russell's December 1992 *Dr. Dobb's Journal* article.
`ALGEBRA.C` is not included in this repository; it is in the article's
companion archive.

```sh
curl -LO http://nemo.asuscomm.com/orbis/geozips/dobbs92.zip
unzip -j dobbs92.zip ALGEBRA.C
mv ALGEBRA.C algebra.c
gcc -O0 -std=c99 -w ref.c -o ref -lm && ./ref
cp stereo.csv ../../tests/data/
```

The pseudo-random sequence is a fixed xorshift, so the output is
reproducible on any IEEE 754 platform with a conforming `libm`.
