      SUBROUTINE CLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                KD, LDAB, N;
      REAL               AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      REAL               S( * );
      COMPLEX            AB( LDAB, * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, THRESH;
      const              ONE = 1.0, THRESH = 0.1 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
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

            // Upper triangle of A is stored in band format.

            for (J = 1; J <= N; J++) { // 20
               CJ = S( J );
               DO 10 I = max( 1, J-KD ), J - 1;
                  AB( KD+1+I-J, J ) = CJ*S( I )*AB( KD+1+I-J, J );
               } // 10
               AB( KD+1, J ) = CJ*CJ*REAL( AB( KD+1, J ) );
            } // 20
         } else {

            // Lower triangle of A is stored.

            for (J = 1; J <= N; J++) { // 40
               CJ = S( J );
               AB( 1, J ) = CJ*CJ*REAL( AB( 1, J ) );
               DO 30 I = J + 1, min( N, J+KD );
                  AB( 1+I-J, J ) = CJ*S( I )*AB( 1+I-J, J );
               } // 30
            } // 40
         }
         EQUED = 'Y';
      }

      return;

      // End of CLAQHB

      }
