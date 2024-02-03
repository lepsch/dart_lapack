      SUBROUTINE CPOTRF( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX            A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      COMPLEX            CONE
      PARAMETER          ( ONE = 1.0E+0, CONE = ( 1.0E+0, 0.0E+0 ) )
*     ..
*     .. Local Scalars ..
      bool               UPPER;
      int                J, JB, NB
*     ..
*     .. External Functions ..
      bool               LSAME;
      int                ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           CGEMM, CHERK, CPOTRF2, CTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPOTRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'CPOTRF', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL CPOTRF2( UPLO, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization A = U**H *U.
*
            DO 10 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL CHERK( 'Upper', 'Conjugate transpose', JB, J-1, -ONE, A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL CPOTRF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block row.
*
                  CALL CGEMM( 'Conjugate transpose', 'No transpose', JB, N-J-JB+1, J-1, -CONE, A( 1, J ), LDA, A( 1, J+JB ), LDA, CONE, A( J, J+JB ), LDA )
                  CALL CTRSM( 'Left', 'Upper', 'Conjugate transpose', 'Non-unit', JB, N-J-JB+1, CONE, A( J, J ), LDA, A( J, J+JB ), LDA )
               END IF
   10       CONTINUE
*
         ELSE
*
*           Compute the Cholesky factorization A = L*L**H.
*
            DO 20 J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL CHERK( 'Lower', 'No transpose', JB, J-1, -ONE, A( J, 1 ), LDA, ONE, A( J, J ), LDA )
               CALL CPOTRF2( 'Lower', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 ) GO TO 30
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block column.
*
                  CALL CGEMM( 'No transpose', 'Conjugate transpose', N-J-JB+1, JB, J-1, -CONE, A( J+JB, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+JB, J ), LDA )
                  CALL CTRSM( 'Right', 'Lower', 'Conjugate transpose', 'Non-unit', N-J-JB+1, JB, CONE, A( J, J ), LDA, A( J+JB, J ), LDA )
               END IF
   20       CONTINUE
         END IF
      END IF
      GO TO 40
*
   30 CONTINUE
      INFO = INFO + J - 1
*
   40 CONTINUE
      RETURN
*
*     End of CPOTRF
*
      END
