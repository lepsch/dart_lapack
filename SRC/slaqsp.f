      SUBROUTINE SLAQSP( UPLO, N, AP, S, SCOND, AMAX, EQUED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                N;
      REAL               AMAX, SCOND
      // ..
      // .. Array Arguments ..
      REAL               AP( * ), S( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, THRESH
      const              ONE = 1.0E+0, THRESH = 0.1E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, JC;
      REAL               CJ, LARGE, SMALL
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF

      // Initialize LARGE and SMALL.

      SMALL = SLAMCH( 'Safe minimum' ) / SLAMCH( 'Precision' )
      LARGE = ONE / SMALL

      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN

         // No equilibration

         EQUED = 'N'
      } else {

         // Replace A by diag(S) * A * diag(S).

         IF( LSAME( UPLO, 'U' ) ) THEN

            // Upper triangle of A is stored.

            JC = 1
            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = 1, J
                  AP( JC+I-1 ) = CJ*S( I )*AP( JC+I-1 )
   10          CONTINUE
               JC = JC + J
   20       CONTINUE
         } else {

            // Lower triangle of A is stored.

            JC = 1
            DO 40 J = 1, N
               CJ = S( J )
               DO 30 I = J, N
                  AP( JC+I-J ) = CJ*S( I )*AP( JC+I-J )
   30          CONTINUE
               JC = JC + N - J + 1
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF

      RETURN

      // End of SLAQSP

      }
