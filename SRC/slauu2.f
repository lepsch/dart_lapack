      SUBROUTINE SLAUU2( UPLO, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      int                INFO, LDA, N
*     ..
*     .. Array Arguments ..
      REAL               A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL               ONE
      PARAMETER          ( ONE = 1.0E+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      int                I
      REAL               AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      REAL               SDOT
      EXTERNAL           LSAME, SDOT
*     ..
*     .. External Subroutines ..
      EXTERNAL           SGEMV, SSCAL, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
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
         CALL XERBLA( 'SLAUU2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Compute the product U * U**T.
*
         DO 10 I = 1, N
            AII = A( I, I )
            IF( I.LT.N ) THEN
               A( I, I ) = SDOT( N-I+1, A( I, I ), LDA, A( I, I ), LDA )
               CALL SGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, AII, A( 1, I ), 1 )
            ELSE
               CALL SSCAL( I, AII, A( 1, I ), 1 )
            END IF
   10    CONTINUE
*
      ELSE
*
*        Compute the product L**T * L.
*
         DO 20 I = 1, N
            AII = A( I, I )
            IF( I.LT.N ) THEN
               A( I, I ) = SDOT( N-I+1, A( I, I ), 1, A( I, I ), 1 )
               CALL SGEMV( 'Transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, AII, A( I, 1 ), LDA )
            ELSE
               CALL SSCAL( I, AII, A( I, 1 ), LDA )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of SLAUU2
*
      END
