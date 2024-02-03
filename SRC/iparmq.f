      int     FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IHI, ILO, ISPEC, LWORK, N;
      String             NAME*( * ), OPTS*( * );

*  ================================================================
      // .. Parameters ..
      int                INMIN, INWIN, INIBL, ISHFTS, IACC22, ICOST;
      const              INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16, ICOST = 17 ;
      int                NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP, RCOST;
      const              NMIN = 75, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500, RCOST = 10 ;
      REAL               TWO
      const              TWO = 2.0 ;
      // ..
      // .. Local Scalars ..
      int                NH, NS;
      int                I, IC, IZ;
      String             SUBNAM*6;
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC LOG, MAX, MOD, NINT, REAL
      // ..
      // .. Executable Statements ..
      if ( ( ISPEC == ISHFTS ) || ( ISPEC == INWIN ) || ( ISPEC == IACC22 ) ) {

         // ==== Set the number simultaneous shifts ====

         NH = IHI - ILO + 1
         NS = 2
         if (NH >= 30) NS = 4          IF( NH >= 60 ) NS = 10          IF( NH >= 150 ) NS = MAX( 10, NH / NINT( LOG( REAL( NH ) ) / LOG( TWO ) ) )          IF( NH >= 590 ) NS = 64          IF( NH >= 3000 ) NS = 128          IF( NH >= 6000 ) NS = 256;
         NS = MAX( 2, NS-MOD( NS, 2 ) )
      }

      if ( ISPEC == INMIN ) {


         // ===== Matrices of order smaller than NMIN get sent
         // .     to xLAHQR, the classic double shift algorithm.
         // .     This must be at least 11. ====

         IPARMQ = NMIN

      } else if ( ISPEC == INIBL ) {

         // ==== INIBL: skip a multi-shift qr iteration and
         // .    whenever aggressive early deflation finds
         // .    at least (NIBBLE*(window size)/100) deflations. ====

         IPARMQ = NIBBLE

      } else if ( ISPEC == ISHFTS ) {

         // ==== NSHFTS: The number of simultaneous shifts =====

         IPARMQ = NS

      } else if ( ISPEC == INWIN ) {

         // ==== NW: deflation window size.  ====

         if ( NH.LE.KNWSWP ) {
            IPARMQ = NS
         } else {
            IPARMQ = 3*NS / 2
         }

      } else if ( ISPEC == IACC22 ) {

         // ==== IACC22: Whether to accumulate reflections
         // .     before updating the far-from-diagonal elements
         // .     and whether to use 2-by-2 block structure while
         // .     doing it.  A small amount of work could be saved
         // .     by making this choice dependent also upon the
         // .     NH=IHI-ILO+1.


         // Convert NAME to upper case if the first character is lower case.

         IPARMQ = 0
         SUBNAM = NAME
         IC = ICHAR( SUBNAM( 1: 1 ) )
         IZ = ICHAR( 'Z' )
         if ( IZ == 90 || IZ == 122 ) {

            // ASCII character set

            if ( IC >= 97 && IC.LE.122 ) {
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               for (I = 2; I <= 6; I++) {
                  IC = ICHAR( SUBNAM( I: I ) )
                  if (IC >= 97 && IC.LE.122) SUBNAM( I: I ) = CHAR( IC-32 );
               }
            }

         } else if ( IZ == 233 || IZ == 169 ) {

            // EBCDIC character set

            if ( ( IC >= 129 && IC.LE.137 ) || ( IC >= 145 && IC.LE.153 ) || ( IC >= 162 && IC.LE.169 ) ) {
               SUBNAM( 1: 1 ) = CHAR( IC+64 )
               for (I = 2; I <= 6; I++) {
                  IC = ICHAR( SUBNAM( I: I ) )
                  IF( ( IC >= 129 && IC.LE.137 ) || ( IC >= 145 && IC.LE.153 ) || ( IC >= 162 && IC.LE.169 ) )SUBNAM( I: I ) = CHAR( IC+64 )
               }
            }

         } else if ( IZ == 218 || IZ == 250 ) {

            // Prime machines:  ASCII+128

            if ( IC >= 225 && IC.LE.250 ) {
               SUBNAM( 1: 1 ) = CHAR( IC-32 )
               for (I = 2; I <= 6; I++) {
                  IC = ICHAR( SUBNAM( I: I ) )
                  if (IC >= 225 && IC.LE.250) SUBNAM( I: I ) = CHAR( IC-32 );
               }
            }
         }

         if ( SUBNAM( 2:6 ) == 'GGHRD' || SUBNAM( 2:6 ) == 'GGHD3' ) {
            IPARMQ = 1
            if (NH >= K22MIN) IPARMQ = 2;
         } else if ( SUBNAM( 4:6 ) == 'EXC' ) {
            if (NH >= KACMIN) IPARMQ = 1             IF( NH >= K22MIN ) IPARMQ = 2          ELSE IF ( SUBNAM( 2:6 ) == 'HSEQR' || SUBNAM( 2:5 ) == 'LAQR' ) THEN             IF( NS >= KACMIN ) IPARMQ = 1             IF( NS >= K22MIN ) IPARMQ = 2;
         }

      } else if ( ISPEC == ICOST ) {

         // === Relative cost of near-the-diagonal chase vs
             // BLAS updates ===

         IPARMQ = RCOST
      } else {
         // ===== invalid value of ispec =====
         IPARMQ = -1

      }

      // ==== End of IPARMQ ====

      }
