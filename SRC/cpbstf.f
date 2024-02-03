      SUBROUTINE CPBSTF( UPLO, N, KD, AB, LDAB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      COMPLEX            AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, KLD, KM, M;
      REAL               AJJ
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHER, CLACGV, CSSCAL, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( KD.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPBSTF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      KLD = MAX( 1, LDAB-1 )

      // Set the splitting point m.

      M = ( N+KD ) / 2

      IF( UPPER ) THEN

         // Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m).

         DO 10 J = N, M + 1, -1

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = REAL( AB( KD+1, J ) )
            IF( AJJ.LE.ZERO ) THEN
               AB( KD+1, J ) = AJJ
               GO TO 50
            END IF
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
            KM = MIN( J-1, KD )

            // Compute elements j-km:j-1 of the j-th column and update the
           t // he leading submatrix within the band.

            CALL CSSCAL( KM, ONE / AJJ, AB( KD+1-KM, J ), 1 )
            CALL CHER( 'Upper', KM, -ONE, AB( KD+1-KM, J ), 1, AB( KD+1, J-KM ), KLD )
   10    CONTINUE

         // Factorize the updated submatrix A(1:m,1:m) as U**H*U.

         DO 20 J = 1, M

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = REAL( AB( KD+1, J ) )
            IF( AJJ.LE.ZERO ) THEN
               AB( KD+1, J ) = AJJ
               GO TO 50
            END IF
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
            KM = MIN( KD, M-J )

            // Compute elements j+1:j+km of the j-th row and update the
           t // railing submatrix within the band.

            IF( KM.GT.0 ) THEN
               CALL CSSCAL( KM, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL CLACGV( KM, AB( KD, J+1 ), KLD )
               CALL CHER( 'Upper', KM, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD )
               CALL CLACGV( KM, AB( KD, J+1 ), KLD )
            END IF
   20    CONTINUE
      ELSE

         // Factorize A(m+1:n,m+1:n) as L**H*L, and update A(1:m,1:m).

         DO 30 J = N, M + 1, -1

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = REAL( AB( 1, J ) )
            IF( AJJ.LE.ZERO ) THEN
               AB( 1, J ) = AJJ
               GO TO 50
            END IF
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
            KM = MIN( J-1, KD )

            // Compute elements j-km:j-1 of the j-th row and update the
           t // railing submatrix within the band.

            CALL CSSCAL( KM, ONE / AJJ, AB( KM+1, J-KM ), KLD )
            CALL CLACGV( KM, AB( KM+1, J-KM ), KLD )
            CALL CHER( 'Lower', KM, -ONE, AB( KM+1, J-KM ), KLD, AB( 1, J-KM ), KLD )
            CALL CLACGV( KM, AB( KM+1, J-KM ), KLD )
   30    CONTINUE

         // Factorize the updated submatrix A(1:m,1:m) as U**H*U.

         DO 40 J = 1, M

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = REAL( AB( 1, J ) )
            IF( AJJ.LE.ZERO ) THEN
               AB( 1, J ) = AJJ
               GO TO 50
            END IF
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
            KM = MIN( KD, M-J )

            // Compute elements j+1:j+km of the j-th column and update the
           t // railing submatrix within the band.

            IF( KM.GT.0 ) THEN
               CALL CSSCAL( KM, ONE / AJJ, AB( 2, J ), 1 )
               CALL CHER( 'Lower', KM, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD )
            END IF
   40    CONTINUE
      END IF
      RETURN

   50 CONTINUE
      INFO = J
      RETURN

      // End of CPBSTF

      }
