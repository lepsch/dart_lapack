      int ilaenv(ISPEC, NAME, OPTS, N1, N2, N3, N4 ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      List<String>       NAME, OPTS;
      int                ISPEC, N1, N2, N3, N4;
      // ..

// =====================================================================

      // .. Intrinsic Functions ..
      // INTRINSIC INT, MIN, REAL
      // ..
      // .. External Functions ..
      int                IEEECK;
      // EXTERNAL IEEECK
      // ..
      // .. Arrays in Common ..
      int                IPARMS( 100 );
      // ..
      // .. Common blocks ..
      // COMMON / CLAENV / IPARMS
      // ..
      // .. Save statement ..
      SAVE               / CLAENV /;
      // ..
      // .. Executable Statements ..

      if ( ISPEC >= 1 && ISPEC <= 5 ) {

         // Return a value from the common block.

         if ( NAME(2:6) == 'GEQR ' ) {
            if (N3 == 2) {
               ILAENV = IPARMS ( 2 );
            } else {
               ILAENV = IPARMS ( 1 );
            }
         } else if ( NAME(2:6) == 'GELQ ' ) {
            if (N3 == 2) {
               ILAENV = IPARMS ( 2 );
            } else {
               ILAENV = IPARMS ( 1 );
            }
         } else {
            ILAENV = IPARMS( ISPEC );
         }

      } else if ( ISPEC == 6 ) {

         // Compute SVD crossover point.

         ILAENV = INT( REAL( min( N1, N2 ) )*1.6 );

      } else if ( ISPEC >= 7 && ISPEC <= 9 ) {

         // Return a value from the common block.

         ILAENV = IPARMS( ISPEC );

      } else if ( ISPEC == 10 ) {

         // IEEE NaN arithmetic can be trusted not to trap

         // ILAENV = 0
         ILAENV = 1;
         if ( ILAENV == 1 ) {
            ILAENV = IEEECK( 1, 0.0, 1.0 );
         }

      } else if ( ISPEC == 11 ) {

         // Infinity arithmetic can be trusted not to trap

         // ILAENV = 0
         ILAENV = 1;
         if ( ILAENV == 1 ) {
            ILAENV = IEEECK( 0, 0.0, 1.0 );
         }

      } else {

         // Invalid value for ISPEC

         ILAENV = -1;
      }

      return;

      // End of ILAENV

      }
      int ilaenv2stage(ISPEC, NAME, OPTS, N1, N2, N3, N4 ) {
      // .. Scalar Arguments ..
      List<String>       NAME, OPTS;
      int                ISPEC, N1, N2, N3, N4;
      // ..

// =====================================================================

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
      SAVE               / CLAENV /;
      // ..
      // .. Executable Statements ..

      if (( ISPEC >= 1 ) && (ISPEC <= 5)) {

      // 1 <= ISPEC <= 5: 2stage eigenvalues SVD routines.

         if ( ISPEC == 1 ) {
             ILAENV2STAGE = IPARMS( 1 );
         } else {
             IISPEC = 16 + ISPEC;
             ILAENV2STAGE = IPARAM2STAGE( IISPEC, NAME, OPTS, N1, N2, N3, N4 );
         }

      } else {

         // Invalid value for ISPEC

         ILAENV2STAGE = -1;
      }

      return;
      }
