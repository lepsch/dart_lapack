      SUBROUTINE CPBTF2( UPLO, N, KD, AB, LDAB, INFO )

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
      PARAMETER          ( ONE = 1.0E+0, ZERO = 0.0E+0 )
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                J, KLD, KN;
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
         CALL XERBLA( 'CPBTF2', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      KLD = MAX( 1, LDAB-1 )

      IF( UPPER ) THEN

         // Compute the Cholesky factorization A = U**H * U.

         DO 10 J = 1, N

            // Compute U(J,J) and test for non-positive-definiteness.

            AJJ = REAL( AB( KD+1, J ) )
            IF( AJJ.LE.ZERO ) THEN
               AB( KD+1, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AB( KD+1, J ) = AJJ

            // Compute elements J+1:J+KN of row J and update the
           t // railing submatrix within the band.

            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL CSSCAL( KN, ONE / AJJ, AB( KD, J+1 ), KLD )
               CALL CLACGV( KN, AB( KD, J+1 ), KLD )
               CALL CHER( 'Upper', KN, -ONE, AB( KD, J+1 ), KLD, AB( KD+1, J+1 ), KLD )
               CALL CLACGV( KN, AB( KD, J+1 ), KLD )
            END IF
   10    CONTINUE
      ELSE

         // Compute the Cholesky factorization A = L*L**H.

         DO 20 J = 1, N

            // Compute L(J,J) and test for non-positive-definiteness.

            AJJ = REAL( AB( 1, J ) )
            IF( AJJ.LE.ZERO ) THEN
               AB( 1, J ) = AJJ
               GO TO 30
            END IF
            AJJ = SQRT( AJJ )
            AB( 1, J ) = AJJ

            // Compute elements J+1:J+KN of column J and update the
           t // railing submatrix within the band.

            KN = MIN( KD, N-J )
            IF( KN.GT.0 ) THEN
               CALL CSSCAL( KN, ONE / AJJ, AB( 2, J ), 1 )
               CALL CHER( 'Lower', KN, -ONE, AB( 2, J ), 1, AB( 1, J+1 ), KLD )
            END IF
   20    CONTINUE
      END IF
      RETURN

   30 CONTINUE
      INFO = J
      RETURN

      // End of CPBTF2

      END
