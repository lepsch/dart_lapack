      SUBROUTINE ZLAQHB( UPLO, N, KD, AB, LDAB, S, SCOND, AMAX, EQUED )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             EQUED, UPLO;
      int                KD, LDAB, N;
      double             AMAX, SCOND;
      // ..
      // .. Array Arguments ..
      double             S( * );
      COMPLEX*16         AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, THRESH;
      const              ONE = 1.0D+0, THRESH = 0.1D+0 ;
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
      // INTRINSIC DBLE, MAX, MIN
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

            // Upper triangle of A is stored in band format.

            DO 20 J = 1, N
               CJ = S( J )
               DO 10 I = MAX( 1, J-KD ), J - 1
                  AB( KD+1+I-J, J ) = CJ*S( I )*AB( KD+1+I-J, J )
   10          CONTINUE
               AB( KD+1, J ) = CJ*CJ*DBLE( AB( KD+1, J ) )
   20       CONTINUE
         ELSE

            // Lower triangle of A is stored.

            DO 40 J = 1, N
               CJ = S( J )
               AB( 1, J ) = CJ*CJ*DBLE( AB( 1, J ) )
               DO 30 I = J + 1, MIN( N, J+KD )
                  AB( 1+I-J, J ) = CJ*S( I )*AB( 1+I-J, J )
   30          CONTINUE
   40       CONTINUE
         END IF
         EQUED = 'Y'
      END IF

      RETURN

      // End of ZLAQHB

      }
