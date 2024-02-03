      int              FUNCTION ILAENV( ISPEC, NAME, OPTS, N1, N2, N3, N4 );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       NAME, OPTS;
      int                ISPEC, N1, N2, N3, N4;
      // ..

*  =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC INT, MIN, REAL
      // ..
      // .. External Functions ..
      int                IEEECK, IPARAM2STAGE;
      // EXTERNAL IEEECK, IPARAM2STAGE
      // ..
      // .. Arrays in Common ..
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / IPARMS
      // ..
      // .. Save statement ..
      SAVE               / CLAENV /
      // ..
      // .. Executable Statements ..

      if ( ISPEC >= 1 && ISPEC.LE.5 ) {

         // Return a value from the common block.

         ILAENV = IPARMS( ISPEC )

      } else if ( ISPEC == 6 ) {

         // Compute SVD crossover point.

         ILAENV = INT( REAL( MIN( N1, N2 ) )*1.6E0 )

      } else if ( ISPEC >= 7 && ISPEC.LE.9 ) {

         // Return a value from the common block.

         ILAENV = IPARMS( ISPEC )

      } else if ( ISPEC == 10 ) {

         // IEEE NaN arithmetic can be trusted not to trap

         // ILAENV = 0
         ILAENV = 1
         if ( ILAENV == 1 ) {
            ILAENV = IEEECK( 1, 0.0, 1.0 )
         }

      } else if ( ISPEC == 11 ) {

         // Infinity arithmetic can be trusted not to trap

         // ILAENV = 0
         ILAENV = 1
         if ( ILAENV == 1 ) {
            ILAENV = IEEECK( 0, 0.0, 1.0 )
         }

      } else if (( ISPEC >= 12 ) && (ISPEC.LE.16)) {

      // 12 <= ISPEC <= 16: xHSEQR or one of its subroutines.

         ILAENV = IPARMS( ISPEC )
          // WRITE(*,*) 'ISPEC = ',ISPEC,' ILAENV =',ILAENV
          // ILAENV = IPARMQ( ISPEC, NAME, OPTS, N1, N2, N3, N4 )

      } else if (( ISPEC >= 17 ) && (ISPEC.LE.21)) {

      // 17 <= ISPEC <= 21: 2stage eigenvalues SVD routines.

         if ( ISPEC == 17 ) {
             ILAENV = IPARMS( 1 )
         } else {
             ILAENV = IPARAM2STAGE( ISPEC, NAME, OPTS, N1, N2, N3, N4 )
         }

      } else {

         // Invalid value for ISPEC

         ILAENV = -1
      }

      RETURN

      // End of ILAENV

      }
      int     FUNCTION ILAENV2STAGE( ISPEC, NAME, OPTS, N1, N2, N3, N4 );
      // .. Scalar Arguments ..
      List<String>       NAME, OPTS;
      int                ISPEC, N1, N2, N3, N4;
      // ..

*  =====================================================================

      // .. Local variables ..
      int                IISPEC;
      // .. External Functions ..
      int                IPARAM2STAGE;
      // EXTERNAL IPARAM2STAGE
      // ..
      // .. Arrays in Common ..
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / IPARMS
      // ..
      // .. Save statement ..
      SAVE               / CLAENV /
      // ..
      // .. Executable Statements ..

      if (( ISPEC >= 1 ) && (ISPEC.LE.5)) {

      // 1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.

         if ( ISPEC == 1 ) {
             ILAENV2STAGE = IPARMS( 1 )
         } else {
             IISPEC = 16 + ISPEC
             ILAENV2STAGE = IPARAM2STAGE( IISPEC, NAME, OPTS, N1, N2, N3, N4 )
         }

      } else {

         // Invalid value for ISPEC

         ILAENV2STAGE = -1
      }

      RETURN

      // End of ILAENV2STAGE

      }
      int     FUNCTION IPARMQ( ISPEC, NAME, OPTS, N, ILO, IHI, LWORK );

      int                INMIN, INWIN, INIBL, ISHFTS, IACC22;
      const              INMIN = 12, INWIN = 13, INIBL = 14, ISHFTS = 15, IACC22 = 16 ;
      int                NMIN, K22MIN, KACMIN, NIBBLE, KNWSWP;
      const              NMIN = 11, K22MIN = 14, KACMIN = 14, NIBBLE = 14, KNWSWP = 500 ;
      REAL               TWO
      const              TWO = 2.0 ;
      // ..
      // .. Scalar Arguments ..
      int                IHI, ILO, ISPEC, LWORK, N;
      String             NAME*( * ), OPTS*( * );
      // ..
      // .. Local Scalars ..
      int                NH, NS;
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
         // .     to LAHQR, the classic double shift algorithm.
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

         IPARMQ = 0
         if (NS >= KACMIN) IPARMQ = 1          IF( NS >= K22MIN ) IPARMQ = 2;

      } else {
         // ===== invalid value of ispec =====
         IPARMQ = -1

      }

      // ==== End of IPARMQ ====

      }
