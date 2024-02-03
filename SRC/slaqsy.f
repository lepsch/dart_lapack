      SUBROUTINE SLAQSY( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                LDA, N;
      REAL               AMAX, SCOND
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), S( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, THRESH
      const              ONE = 1.0E+0, THRESH = 0.1E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               CJ, LARGE, SMALL
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.0 ) {
         EQUED = 'N'
         RETURN
      }

      // Initialize LARGE and SMALL.

      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      LARGE = ONE / SMALL

      if ( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) {

         // No equilibration

         EQUED = 'N'
      } else {

         // Replace A by diag(S) * A * diag(S).

         if ( LSAME( UPLO, 'U' ) ) {

            // Upper triangle of A is stored.

            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = 1, J
                  A( I, J ) = CJ*S( I )*A( I, J )
   10          CONTINUE
   20       CONTINUE
         } else {

            // Lower triangle of A is stored.

            DO 40 J = 1, N
               CJ = S( J )
               DO 30 I = J, N
                  A( I, J ) = CJ*S( I )*A( I, J )
   30          CONTINUE
   40       CONTINUE
         }
         EQUED = 'Y'
      }

      RETURN

      // End of SLAQSY

      }
