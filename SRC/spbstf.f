      SUBROUTINE SPBSTF( UPLO, N, KD, AB, LDAB, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, KD, LDAB, N;
      // ..
      // .. Array Arguments ..
      REAL               AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
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
      // EXTERNAL SSCAL, SSYR, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, SQRT
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
         CALL XERBLA( 'SPBSTF', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      KLD = MAX( 1, LDAB-1 )

      // Set the splitting point m.

      M = ( N+KD ) / 2

      IF( UPPER ) THEN

         // Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).

         DO 10 J = N, M + 1, -1

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
            KM = MIN( J-1, KD )

            // Compute elements j-km:j-1 of the j-th column and update the
           t // he leading submatrix within the band.

            CALL SSCAL( KM, ONE / AJJ, AB( KD+1-KM, J ), 1 )
            CALL SSYR( 'Upper', KM, -ONE, AB( KD+1-KM, J ), 1, AB( KD+1, J-KM ), KLD )
   10    CONTINUE

         // Factorize the updated submatrix A(1:m,1:m) as U**T*U.

         DO 20 J = 1, M

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( KD+1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ
            KM = MIN( KD, M-J )

            // Compute elements j+1:j+km of the j-th row and update the
           t // railing submatrix within the band.

            IF( KM.GT.0 ) THEN
               CALL SSCAL( KM, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL SSYR( 'Upper', KM, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD )
            END IF
   20    CONTINUE
      ELSE

         // Factorize A(m+1:n,m+1:n) as L**T*L, and update A(1:m,1:m).

         DO 30 J = N, M + 1, -1

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
            KM = MIN( J-1, KD )

            // Compute elements j-km:j-1 of the j-th row and update the
           t // railing submatrix within the band.

            CALL SSCAL( KM, ONE / AJJ, AB( KM+1, J-KM ), KLD )
            CALL SSYR( 'Lower', KM, -ONE, AB( KM+1, J-KM ), KLD, AB( 1, J-KM ), KLD )
   30    CONTINUE

         // Factorize the updated submatrix A(1:m,1:m) as U**T*U.

         DO 40 J = 1, M

            // Compute s(j,j) and test for non-positive-definiteness.

            AJJ = AB( 1, J )
            IF( AJJ.LE.ZERO ) GO TO 50
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ
            KM = MIN( KD, M-J )

            // Compute elements j+1:j+km of the j-th column and update the
           t // railing submatrix within the band.

            IF( KM.GT.0 ) THEN
               CALL SSCAL( KM, ONE / AJJ, AB( 2, J ), 1 )
               CALL SSYR( 'Lower', KM, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD )
            END IF
   40    CONTINUE
      END IF
      RETURN

   50 CONTINUE
      INFO = J
      RETURN

      // End of SPBSTF

      END
