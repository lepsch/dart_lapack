      SUBROUTINE CLAQHP( UPLO, N, AP, S, SCOND, AMAX, EQUED );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                N;
      REAL               AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      REAL               S( * );
      COMPLEX            AP( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JC;
      REAL               CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC REAL
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N <= 0 ) {
         EQUED = 'N';
         return;
      }

      // Initialize LARGE and SMALL.

      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' );
      LARGE = ONE / SMALL;

      if ( SCOND >= THRESH && AMAX >= SMALL && AMAX <= LARGE ) {

         // No equilibration

         EQUED = 'N';
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( LSAME( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            JC = 1;
            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               for (I = 1; I <= J - 1; I++) { // 10
                  AP( JC+I-1 ) = CJ*S( I )*AP( JC+I-1 );
               } // 10
               AP( JC+J-1 ) = CJ*CJ*REAL( AP( JC+J-1 ) );
               JC = JC + J;
            } // 20
         } else {

            // Lower triangle of A is stored.

            JC = 1;
            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               AP( JC ) = CJ*CJ*REAL( AP( JC ) );
               for (I = J + 1; I <= N; I++) { // 30
                  AP( JC+I-J ) = CJ*S( I )*AP( JC+I-J );
               } // 30
               JC = JC + N - J + 1;
            } // 40
         }
         EQUED = 'Y';
      }

      return;

      // End of CLAQHP

      }
