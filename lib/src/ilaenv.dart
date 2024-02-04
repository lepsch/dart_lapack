import 'dart:math';

import 'package:lapack/src/ieeeck.dart';
import 'package:lapack/src/iparmq.dart';

int ilaenv(
  final int ISPEC,
  final String NAME,
  final String OPTS,
  final int N1,
  final int N2,
  final int N3,
  final int N4,
) {
// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
  int NB, NBMIN, NX;
  bool CNAME, SNAME, TWOSTAGE;
  String C1 = '', C2 = '', C4 = '', C3 = '', SUBNAM = '';

  switch (ISPEC) {
    case 1:
    case 2:
    case 3:
      // Convert NAME to upper case in case the first character is lower case.
      SUBNAM = NAME.toUpperCase();

      C1 = SUBNAM.substring(0, 1);
      SNAME = C1 == 'S' || C1 == 'D';
      CNAME = C1 == 'C' || C1 == 'Z';
      if (!(CNAME || SNAME)) return 1;
      C2 = SUBNAM.substring(1, 3);
      C3 = SUBNAM.substring(3, 6);
      C4 = C3.substring(1, 3);
      TWOSTAGE = SUBNAM.length >= 11 && SUBNAM.substring(10, 11) == '2';

      if (ISPEC == 1) {
        // ISPEC = 1:  block size

        // In these examples, separate code is provided for setting NB for
        // real and complex.  We assume that NB will take the same value in
        // single or double precision.

        NB = 1;

        if (SUBNAM.substring(1, 6) == 'LAORH') {
          // This is for *LAORHR_GETRFNP routine

          if (SNAME) {
            NB = 32;
          } else {
            NB = 32;
          }
        } else if (C2 == 'GE') {
          if (C3 == 'TRF') {
            if (SNAME) {
              NB = 64;
            } else {
              NB = 64;
            }
          } else if (C3 == 'QRF' || C3 == 'RQF' || C3 == 'LQF' || C3 == 'QLF') {
            if (SNAME) {
              NB = 32;
            } else {
              NB = 32;
            }
          } else if (C3 == 'QR ') {
            if (N3 == 1) {
              if (SNAME) {
                // M*N
                if ((N1 * N2 <= 131072) || (N1 <= 8192)) {
                  NB = N1;
                } else {
                  NB = 32768 ~/ N2;
                }
              } else {
                if ((N1 * N2 <= 131072) || (N1 <= 8192)) {
                  NB = N1;
                } else {
                  NB = 32768 ~/ N2;
                }
              }
            } else {
              if (SNAME) {
                NB = 1;
              } else {
                NB = 1;
              }
            }
          } else if (C3 == 'LQ ') {
            if (N3 == 2) {
              if (SNAME) {
                // M*N
                if ((N1 * N2 <= 131072) || (N1 <= 8192)) {
                  NB = N1;
                } else {
                  NB = 32768 ~/ N2;
                }
              } else {
                if ((N1 * N2 <= 131072) || (N1 <= 8192)) {
                  NB = N1;
                } else {
                  NB = 32768 ~/ N2;
                }
              }
            } else {
              if (SNAME) {
                NB = 1;
              } else {
                NB = 1;
              }
            }
          } else if (C3 == 'HRD') {
            if (SNAME) {
              NB = 32;
            } else {
              NB = 32;
            }
          } else if (C3 == 'BRD') {
            if (SNAME) {
              NB = 32;
            } else {
              NB = 32;
            }
          } else if (C3 == 'TRI') {
            if (SNAME) {
              NB = 64;
            } else {
              NB = 64;
            }
          } else if (SUBNAM.substring(3, 7) == 'QP3RK') {
            if (SNAME) {
              NB = 32;
            } else {
              NB = 32;
            }
          }
        } else if (C2 == 'PO') {
          if (C3 == 'TRF') {
            if (SNAME) {
              NB = 64;
            } else {
              NB = 64;
            }
          }
        } else if (C2 == 'SY') {
          if (C3 == 'TRF') {
            if (SNAME) {
              if (TWOSTAGE) {
                NB = 192;
              } else {
                NB = 64;
              }
            } else {
              if (TWOSTAGE) {
                NB = 192;
              } else {
                NB = 64;
              }
            }
          } else if (SNAME && C3 == 'TRD') {
            NB = 32;
          } else if (SNAME && C3 == 'GST') {
            NB = 64;
          }
        } else if (CNAME && C2 == 'HE') {
          if (C3 == 'TRF') {
            if (TWOSTAGE) {
              NB = 192;
            } else {
              NB = 64;
            }
          } else if (C3 == 'TRD') {
            NB = 32;
          } else if (C3 == 'GST') {
            NB = 64;
          }
        } else if (SNAME && C2 == 'OR') {
          if (C3.substring(0, 1) == 'G') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NB = 32;
            }
          } else if (C3.substring(0, 1) == 'M') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NB = 32;
            }
          }
        } else if (CNAME && C2 == 'UN') {
          if (C3.substring(0, 1) == 'G') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NB = 32;
            }
          } else if (C3.substring(0, 1) == 'M') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NB = 32;
            }
          }
        } else if (C2 == 'GB') {
          if (C3 == 'TRF') {
            if (SNAME) {
              if (N4 <= 64) {
                NB = 1;
              } else {
                NB = 32;
              }
            } else {
              if (N4 <= 64) {
                NB = 1;
              } else {
                NB = 32;
              }
            }
          }
        } else if (C2 == 'PB') {
          if (C3 == 'TRF') {
            if (SNAME) {
              if (N2 <= 64) {
                NB = 1;
              } else {
                NB = 32;
              }
            } else {
              if (N2 <= 64) {
                NB = 1;
              } else {
                NB = 32;
              }
            }
          }
        } else if (C2 == 'TR') {
          if (C3 == 'TRI') {
            if (SNAME) {
              NB = 64;
            } else {
              NB = 64;
            }
          } else if (C3 == 'EVC') {
            if (SNAME) {
              NB = 64;
            } else {
              NB = 64;
            }
          } else if (C3 == 'SYL') {
            // The upper bound is to prevent overly aggressive scaling.
            if (SNAME) {
              NB = min(max(48, (min(N1, N2) * 16) ~/ 100), 240);
            } else {
              NB = min(max(24, (min(N1, N2) * 8) ~/ 100), 80);
            }
          }
        } else if (C2 == 'LA') {
          if (C3 == 'UUM') {
            if (SNAME) {
              NB = 64;
            } else {
              NB = 64;
            }
          } else if (C3 == 'TRS') {
            if (SNAME) {
              NB = 32;
            } else {
              NB = 32;
            }
          }
        } else if (SNAME && C2 == 'ST') {
          if (C3 == 'EBZ') {
            NB = 1;
          }
        } else if (C2 == 'GG') {
          NB = 32;
          if (C3 == 'HD3') {
            if (SNAME) {
              NB = 32;
            } else {
              NB = 32;
            }
          }
        }
        return NB;
      }

      if (ISPEC == 2) {
        // ISPEC = 2:  minimum block size

        NBMIN = 2;
        if (C2 == 'GE') {
          if (C3 == 'QRF' || C3 == 'RQF' || C3 == 'LQF' || C3 == 'QLF') {
            if (SNAME) {
              NBMIN = 2;
            } else {
              NBMIN = 2;
            }
          } else if (C3 == 'HRD') {
            if (SNAME) {
              NBMIN = 2;
            } else {
              NBMIN = 2;
            }
          } else if (C3 == 'BRD') {
            if (SNAME) {
              NBMIN = 2;
            } else {
              NBMIN = 2;
            }
          } else if (C3 == 'TRI') {
            if (SNAME) {
              NBMIN = 2;
            } else {
              NBMIN = 2;
            }
          } else if (SUBNAM.substring(3, 7) == 'QP3RK') {
            if (SNAME) {
              NBMIN = 2;
            } else {
              NBMIN = 2;
            }
          }
        } else if (C2 == 'SY') {
          if (C3 == 'TRF') {
            if (SNAME) {
              NBMIN = 8;
            } else {
              NBMIN = 8;
            }
          } else if (SNAME && C3 == 'TRD') {
            NBMIN = 2;
          }
        } else if (CNAME && C2 == 'HE') {
          if (C3 == 'TRD') {
            NBMIN = 2;
          }
        } else if (SNAME && C2 == 'OR') {
          if (C3.substring(0, 1) == 'G') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NBMIN = 2;
            }
          } else if (C3.substring(0, 1) == 'M') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NBMIN = 2;
            }
          }
        } else if (CNAME && C2 == 'UN') {
          if (C3.substring(0, 1) == 'G') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NBMIN = 2;
            }
          } else if (C3.substring(0, 1) == 'M') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NBMIN = 2;
            }
          }
        } else if (C2 == 'GG') {
          NBMIN = 2;
          if (C3 == 'HD3') {
            NBMIN = 2;
          }
        }
        return NBMIN;
      }

      if (ISPEC == 3) {
        // ISPEC = 3:  crossover point
        NX = 0;
        if (C2 == 'GE') {
          if (C3 == 'QRF' || C3 == 'RQF' || C3 == 'LQF' || C3 == 'QLF') {
            if (SNAME) {
              NX = 128;
            } else {
              NX = 128;
            }
          } else if (C3 == 'HRD') {
            if (SNAME) {
              NX = 128;
            } else {
              NX = 128;
            }
          } else if (C3 == 'BRD') {
            if (SNAME) {
              NX = 128;
            } else {
              NX = 128;
            }
          } else if (SUBNAM.substring(3, 7) == 'QP3RK') {
            if (SNAME) {
              NX = 128;
            } else {
              NX = 128;
            }
          }
        } else if (C2 == 'SY') {
          if (SNAME && C3 == 'TRD') {
            NX = 32;
          }
        } else if (CNAME && C2 == 'HE') {
          if (C3 == 'TRD') {
            NX = 32;
          }
        } else if (SNAME && C2 == 'OR') {
          if (C3.substring(0, 1) == 'G') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NX = 128;
            }
          }
        } else if (CNAME && C2 == 'UN') {
          if (C3.substring(0, 1) == 'G') {
            if (C4 == 'QR' ||
                C4 == 'RQ' ||
                C4 == 'LQ' ||
                C4 == 'QL' ||
                C4 == 'HR' ||
                C4 == 'TR' ||
                C4 == 'BR') {
              NX = 128;
            }
          }
        } else if (C2 == 'GG') {
          NX = 128;
          if (C3 == 'HD3') {
            NX = 128;
          }
        }
        return NX;
      }

    case 4:
      // ISPEC = 4:  number of shifts (used by xHSEQR)
      return 6;

    case 5:
      // ISPEC = 5:  minimum column dimension (not used)
      return 2;

    case 6:
      // ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)
      return (min(N1, N2) * 1.6).toInt();

    case 7:
      // ISPEC = 7:  number of processors (not used)
      return 1;

    case 8:
      // ISPEC = 8:  crossover point for multishift (used by xHSEQR)
      return 50;

    case 9:
      // ISPEC = 9:  maximum size of the subproblems at the bottom of the
      // computation tree in the divide-and-conquer algorithm
      // (used by xGELSD and xGESDD)
      return 25;

    case 10:
      // ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap
      return ieeeck(1, 0.0, 1.0);

    case 11:

      // ISPEC = 11: ieee infinity arithmetic can be trusted not to trap
      return ieeeck(0, 0.0, 1.0);

    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
      // 12 <= ISPEC <= 17: xHSEQR or related subroutines.
      return iparmq(ISPEC, NAME, OPTS, N1, N2, N3, N4);
  }
  // Invalid value for ISPEC
  return -1;
}
