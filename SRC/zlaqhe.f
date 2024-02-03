      SUBROUTINE ZLAQHE( UPLO, N, A, LDA, S, SCOND, AMAX, EQUED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                LDA, N;
      double             AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      double             S( * );
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, THRESH;
      PARAMETER          ( ONE = 1.0D+0, THRESH = 0.1D+0 )
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             CJ, LARGE, SMALL;
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.LE.0 ) THEN
         EQUED = 'N'
         RETURN
      END IF

      // Initialize LARGE and SMALL.

      SMALL = DLAMCH( 'Safe minimum' ) / DLAMCH( 'Precision' )
      LARGE = ONE / SMALL

      IF( SCOND.GE.THRESH .AND. AMAX.GE.SMALL .AND. AMAX.LE.LARGE ) THEN

         // No equilibration

         EQUED = 'N'
      ELSE

         // Replace A by diag(S) * A * diag(S).

         IF( LSAME( UPLO, 'U' ) ) THEN

            // Upper triangle of A is stored.

            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = 1, J - 1
                  A( I, J ) = CJ*S( I )*A( I, J )
   10          CONTINUE
               A( J, J ) = CJ*CJ*DBLE( A( J, J ) )
   20       CONTINUE
         ELSE

            // Lower triangle of A is stored.

            DO 40 J = 1, N
               CJ = S( J )
               A( J, J ) = CJ*CJ*DBLE( A( J, J ) )
               DO 30 I = J + 1, N
                  A( I, J ) = CJ*S( I )*A( I, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF

      RETURN

      // End of ZLAQHE

      }
