      SUBROUTINE SLARRA( N, D, E, E2, SPLTOL, TNRM, NSPLIT, ISPLIT, INFO );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N, NSPLIT;
      REAL                SPLTOL, TNRM;
      // ..
      // .. Array Arguments ..
      int                ISPLIT( * );
      REAL               D( * ), E( * ), E2( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      REAL               EABS, TMP1;

      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      INFO = 0;
      NSPLIT = 1;

      // Quick return if possible

      if ( N <= 0 ) {
         return;
      }

      // Compute splitting points
      if (SPLTOL < ZERO) {
         // Criterion based on absolute off-diagonal value
         TMP1 = ABS(SPLTOL)* TNRM;
         for (I = 1; I <= N-1; I++) { // 9
            EABS = ABS( E(I) );
            if ( EABS <= TMP1) {
               E(I) = ZERO;
               E2(I) = ZERO;
               ISPLIT( NSPLIT ) = I;
               NSPLIT = NSPLIT + 1;
            }
         } // 9
      } else {
         // Criterion that guarantees relative accuracy
         for (I = 1; I <= N-1; I++) { // 10
            EABS = ABS( E(I) );
            if ( EABS <= SPLTOL * sqrt(ABS(D(I)))*sqrt(ABS(D(I+1))) ) {
               E(I) = ZERO;
               E2(I) = ZERO;
               ISPLIT( NSPLIT ) = I;
               NSPLIT = NSPLIT + 1;
            }
         } // 10
      }
      ISPLIT( NSPLIT ) = N;

      return;

      // End of SLARRA

      }
