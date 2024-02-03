      int ilaenv2stage(ISPEC, NAME, OPTS, N1, N2, N3, N4 ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      // July 2017

      // .. Scalar Arguments ..
      List<String>       NAME, OPTS;
      int                ISPEC, N1, N2, N3, N4;
      // ..

// =====================================================================
      // ..
      // .. Local Scalars ..
      int                IISPEC;
      // ..
      // .. External Functions ..
      //- int                IPARAM2STAGE;
      // EXTERNAL IPARAM2STAGE
      // ..
      // .. Executable Statements ..

      GO TO ( 10, 10, 10, 10, 10 )ISPEC;

      // Invalid value for ISPEC

      ILAENV2STAGE = -1;
      return;

      } // 10

      // 2stage eigenvalues and SVD or related subroutines.

      IISPEC = 16 + ISPEC;
      ILAENV2STAGE = IPARAM2STAGE( IISPEC, NAME, OPTS, N1, N2, N3, N4 );
      return;
      }
