      int     FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       NAME, OPTS;
      int                ISPEC, N1, N2, N3, N4;
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IC, IZ, NB, NBMIN, NX;
      bool               CNAME, SNAME, TWOSTAGE;
      String             C1*1, C2*2, C4*2, C3*3, SUBNAM*16;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CHAR, ICHAR, INT, MIN, REAL
      // ..
      // .. External Functions ..
      int                IEEECK, IPARMQ, IPARAM2STAGE;
      // EXTERNAL IEEECK, IPARMQ, IPARAM2STAGE
      // ..
      // .. Executable Statements ..

      GO TO ( 10, 10, 10, 80, 90, 100, 110, 120, 130, 140, 150, 160, 160, 160, 160, 160, 160)ISPEC;

      // Invalid value for ISPEC

      ILAENV = -1;
      return;

      } // 10

      // Convert NAME to upper case if the first character is lower case.

      ILAENV = 1;
      SUBNAM = NAME;
      IC = ICHAR( SUBNAM( 1: 1 ) );
      IZ = ICHAR( 'Z' );
      if ( IZ == 90 || IZ == 122 ) {

         // ASCII character set

         if ( IC >= 97 && IC <= 122 ) {
            SUBNAM( 1: 1 ) = CHAR( IC-32 );
            for (I = 2; I <= 6; I++) { // 20
               IC = ICHAR( SUBNAM( I: I ) );
               if (IC >= 97 && IC <= 122) SUBNAM( I: I ) = CHAR( IC-32 );
            } // 20
         }

      } else if ( IZ == 233 || IZ == 169 ) {

         // EBCDIC character set

         if ( ( IC >= 129 && IC <= 137 ) || ( IC >= 145 && IC <= 153 ) || ( IC >= 162 && IC <= 169 ) ) {
            SUBNAM( 1: 1 ) = CHAR( IC+64 );
            for (I = 2; I <= 6; I++) { // 30
               IC = ICHAR( SUBNAM( I: I ) );
               if( ( IC >= 129 && IC <= 137 ) || ( IC >= 145 && IC <= 153 ) || ( IC >= 162 && IC <= 169 ) )SUBNAM( I: I ) = CHAR( IC+64 );
            } // 30
         }

      } else if ( IZ == 218 || IZ == 250 ) {

         // Prime machines:  ASCII+128

         if ( IC >= 225 && IC <= 250 ) {
            SUBNAM( 1: 1 ) = CHAR( IC-32 );
            for (I = 2; I <= 6; I++) { // 40
               IC = ICHAR( SUBNAM( I: I ) );
               if (IC >= 225 && IC <= 250) SUBNAM( I: I ) = CHAR( IC-32 );
            } // 40
         }
      }

      C1 = SUBNAM( 1: 1 );
      SNAME = C1 == 'S' || C1 == 'D';
      CNAME = C1 == 'C' || C1 == 'Z';
      if( !( CNAME || SNAME ) ) return;
      C2 = SUBNAM( 2: 3 );
      C3 = SUBNAM( 4: 6 );
      C4 = C3( 2: 3 );
      TWOSTAGE = LEN( SUBNAM ) >= 11 && SUBNAM( 11: 11 ) == '2';

      GO TO ( 50, 60, 70 )ISPEC;

      } // 50

      // ISPEC = 1:  block size

      // In these examples, separate code is provided for setting NB for
      // real and complex.  We assume that NB will take the same value in
      // single or double precision.

      NB = 1;

      if ( SUBNAM(2:6) == 'LAORH' ) {

         // This is for *LAORHR_GETRFNP routine

         if ( SNAME ) {
             NB = 32;
         } else {
             NB = 32;
         }
      } else if ( C2 == 'GE' ) {
         if ( C3 == 'TRF' ) {
            if ( SNAME ) {
               NB = 64;
            } else {
               NB = 64;
            }
         } else if ( C3 == 'QRF' || C3 == 'RQF' || C3 == 'LQF' || C3 == 'QLF' ) {
            if ( SNAME ) {
               NB = 32;
            } else {
               NB = 32;
            }
         } else if ( C3 == 'QR ') {
            if ( N3 == 1) {
               if ( SNAME ) {
      // M*N
                  if ((N1*N2 <= 131072) || (N1 <= 8192)) {
                     NB = N1;
                  } else {
                     NB = 32768/N2;
                  }
               } else {
                  if ((N1*N2 <= 131072) || (N1 <= 8192)) {
                     NB = N1;
                  } else {
                     NB = 32768/N2;
                  }
               }
            } else {
               if ( SNAME ) {
                  NB = 1;
               } else {
                  NB = 1;
               }
            }
         } else if ( C3 == 'LQ ') {
            if ( N3 == 2) {
               if ( SNAME ) {
      // M*N
                  if ((N1*N2 <= 131072) || (N1 <= 8192)) {
                     NB = N1;
                  } else {
                     NB = 32768/N2;
                  }
               } else {
                  if ((N1*N2 <= 131072) || (N1 <= 8192)) {
                     NB = N1;
                  } else {
                     NB = 32768/N2;
                  }
               }
            } else {
               if ( SNAME ) {
                  NB = 1;
               } else {
                  NB = 1;
               }
            }
         } else if ( C3 == 'HRD' ) {
            if ( SNAME ) {
               NB = 32;
            } else {
               NB = 32;
            }
         } else if ( C3 == 'BRD' ) {
            if ( SNAME ) {
               NB = 32;
            } else {
               NB = 32;
            }
         } else if ( C3 == 'TRI' ) {
            if ( SNAME ) {
               NB = 64;
            } else {
               NB = 64;
            }
         } else if ( SUBNAM( 4: 7 ) == 'QP3RK' ) {
            if ( SNAME ) {
               NB = 32;
            } else {
               NB = 32;
            }
         }
      } else if ( C2 == 'PO' ) {
         if ( C3 == 'TRF' ) {
            if ( SNAME ) {
               NB = 64;
            } else {
               NB = 64;
            }
         }
      } else if ( C2 == 'SY' ) {
         if ( C3 == 'TRF' ) {
            if ( SNAME ) {
               if ( TWOSTAGE ) {
                  NB = 192;
               } else {
                  NB = 64;
               }
            } else {
               if ( TWOSTAGE ) {
                  NB = 192;
               } else {
                  NB = 64;
               }
            }
         } else if ( SNAME && C3 == 'TRD' ) {
            NB = 32;
         } else if ( SNAME && C3 == 'GST' ) {
            NB = 64;
         }
      } else if ( CNAME && C2 == 'HE' ) {
         if ( C3 == 'TRF' ) {
            if ( TWOSTAGE ) {
               NB = 192;
            } else {
               NB = 64;
            }
         } else if ( C3 == 'TRD' ) {
            NB = 32;
         } else if ( C3 == 'GST' ) {
            NB = 64;
         }
      } else if ( SNAME && C2 == 'OR' ) {
         if ( C3( 1: 1 ) == 'G' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NB = 32;
            }
         } else if ( C3( 1: 1 ) == 'M' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NB = 32;
            }
         }
      } else if ( CNAME && C2 == 'UN' ) {
         if ( C3( 1: 1 ) == 'G' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NB = 32;
            }
         } else if ( C3( 1: 1 ) == 'M' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NB = 32;
            }
         }
      } else if ( C2 == 'GB' ) {
         if ( C3 == 'TRF' ) {
            if ( SNAME ) {
               if ( N4 <= 64 ) {
                  NB = 1;
               } else {
                  NB = 32;
               }
            } else {
               if ( N4 <= 64 ) {
                  NB = 1;
               } else {
                  NB = 32;
               }
            }
         }
      } else if ( C2 == 'PB' ) {
         if ( C3 == 'TRF' ) {
            if ( SNAME ) {
               if ( N2 <= 64 ) {
                  NB = 1;
               } else {
                  NB = 32;
               }
            } else {
               if ( N2 <= 64 ) {
                  NB = 1;
               } else {
                  NB = 32;
               }
            }
         }
      } else if ( C2 == 'TR' ) {
         if ( C3 == 'TRI' ) {
            if ( SNAME ) {
               NB = 64;
            } else {
               NB = 64;
            }
         } else if ( C3 == 'EVC' ) {
            if ( SNAME ) {
               NB = 64;
            } else {
               NB = 64;
            }
         } else if ( C3 == 'SYL' ) {
            // The upper bound is to prevent overly aggressive scaling.
            if ( SNAME ) {
               NB = MIN( MAX( 48, INT( ( MIN( N1, N2 ) * 16 ) / 100) ), 240 );
            } else {
               NB = MIN( MAX( 24, INT( ( MIN( N1, N2 ) * 8 ) / 100) ), 80 );
            }
         }
      } else if ( C2 == 'LA' ) {
         if ( C3 == 'UUM' ) {
            if ( SNAME ) {
               NB = 64;
            } else {
               NB = 64;
            }
         } else if ( C3 == 'TRS' ) {
            if ( SNAME ) {
               NB = 32;
            } else {
               NB = 32;
            }
         }
      } else if ( SNAME && C2 == 'ST' ) {
         if ( C3 == 'EBZ' ) {
            NB = 1;
         }
      } else if ( C2 == 'GG' ) {
         NB = 32;
         if ( C3 == 'HD3' ) {
            if ( SNAME ) {
               NB = 32;
            } else {
               NB = 32;
            }
         }
      }
      ILAENV = NB;
      return;

      } // 60

      // ISPEC = 2:  minimum block size

      NBMIN = 2;
      if ( C2 == 'GE' ) {
         if ( C3 == 'QRF' || C3 == 'RQF' || C3 == 'LQF' || C3 == 'QLF' ) {
            if ( SNAME ) {
               NBMIN = 2;
            } else {
               NBMIN = 2;
            }
         } else if ( C3 == 'HRD' ) {
            if ( SNAME ) {
               NBMIN = 2;
            } else {
               NBMIN = 2;
            }
         } else if ( C3 == 'BRD' ) {
            if ( SNAME ) {
               NBMIN = 2;
            } else {
               NBMIN = 2;
            }
         } else if ( C3 == 'TRI' ) {
            if ( SNAME ) {
               NBMIN = 2;
            } else {
               NBMIN = 2;
            }
         } else if ( SUBNAM( 4: 7 ) == 'QP3RK' ) {
            if ( SNAME ) {
               NBMIN = 2;
            } else {
               NBMIN = 2;
            }
         }

      } else if ( C2 == 'SY' ) {
         if ( C3 == 'TRF' ) {
            if ( SNAME ) {
               NBMIN = 8;
            } else {
               NBMIN = 8;
            }
         } else if ( SNAME && C3 == 'TRD' ) {
            NBMIN = 2;
         }
      } else if ( CNAME && C2 == 'HE' ) {
         if ( C3 == 'TRD' ) {
            NBMIN = 2;
         }
      } else if ( SNAME && C2 == 'OR' ) {
         if ( C3( 1: 1 ) == 'G' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NBMIN = 2;
            }
         } else if ( C3( 1: 1 ) == 'M' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NBMIN = 2;
            }
         }
      } else if ( CNAME && C2 == 'UN' ) {
         if ( C3( 1: 1 ) == 'G' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NBMIN = 2;
            }
         } else if ( C3( 1: 1 ) == 'M' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NBMIN = 2;
            }
         }
      } else if ( C2 == 'GG' ) {
         NBMIN = 2;
         if ( C3 == 'HD3' ) {
            NBMIN = 2;
         }
      }
      ILAENV = NBMIN;
      return;

      } // 70

      // ISPEC = 3:  crossover point

      NX = 0;
      if ( C2 == 'GE' ) {
         if ( C3 == 'QRF' || C3 == 'RQF' || C3 == 'LQF' || C3 == 'QLF' ) {
            if ( SNAME ) {
               NX = 128;
            } else {
               NX = 128;
            }
         } else if ( C3 == 'HRD' ) {
            if ( SNAME ) {
               NX = 128;
            } else {
               NX = 128;
            }
         } else if ( C3 == 'BRD' ) {
            if ( SNAME ) {
               NX = 128;
            } else {
               NX = 128;
            }
         } else if ( SUBNAM( 4: 7 ) == 'QP3RK' ) {
            if ( SNAME ) {
               NX = 128;
            } else {
               NX = 128;
            }
         }
      } else if ( C2 == 'SY' ) {
         if ( SNAME && C3 == 'TRD' ) {
            NX = 32;
         }
      } else if ( CNAME && C2 == 'HE' ) {
         if ( C3 == 'TRD' ) {
            NX = 32;
         }
      } else if ( SNAME && C2 == 'OR' ) {
         if ( C3( 1: 1 ) == 'G' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NX = 128;
            }
         }
      } else if ( CNAME && C2 == 'UN' ) {
         if ( C3( 1: 1 ) == 'G' ) {
            if ( C4 == 'QR' || C4 == 'RQ' || C4 == 'LQ' || C4 == 'QL' || C4 == 'HR' || C4 == 'TR' || C4 == 'BR' ) {
               NX = 128;
            }
         }
      } else if ( C2 == 'GG' ) {
         NX = 128;
         if ( C3 == 'HD3' ) {
            NX = 128;
         }
      }
      ILAENV = NX;
      return;

      } // 80

      // ISPEC = 4:  number of shifts (used by xHSEQR)

      ILAENV = 6;
      return;

      } // 90

      // ISPEC = 5:  minimum column dimension (not used)

      ILAENV = 2;
      return;

      } // 100

      // ISPEC = 6:  crossover point for SVD (used by xGELSS and xGESVD)

      ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6 );
      return;

      } // 110

      // ISPEC = 7:  number of processors (not used)

      ILAENV = 1;
      return;

      } // 120

      // ISPEC = 8:  crossover point for multishift (used by xHSEQR)

      ILAENV = 50;
      return;

      } // 130

      // ISPEC = 9:  maximum size of the subproblems at the bottom of the
                  // computation tree in the divide-and-conquer algorithm
                  // (used by xGELSD and xGESDD)

      ILAENV = 25;
      return;

      } // 140

      // ISPEC = 10: ieee and infinity NaN arithmetic can be trusted not to trap

      // ILAENV = 0
      ILAENV = 1;
      if ( ILAENV == 1 ) {
         ILAENV = IEEECK( 1, 0.0, 1.0 );
      }
      return;

      } // 150

      // ISPEC = 11: ieee infinity arithmetic can be trusted not to trap

      // ILAENV = 0
      ILAENV = 1;
      if ( ILAENV == 1 ) {
         ILAENV = IEEECK( 0, 0.0, 1.0 );
      }
      return;

      } // 160

      // 12 <= ISPEC <= 17: xHSEQR or related subroutines.

      ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 );
      return;

      // End of ILAENV

      }
