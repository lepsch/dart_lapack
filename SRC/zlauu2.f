      SUBROUTINE ZLAUU2( UPLO, N, A, LDA, INFO )
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
      COMPLEX*16         A( LDA, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      int                I
      DOUBLE PRECISION   AII
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZGEMV, ZLACGV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, MAX
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
         CALL XERBLA( 'ZLAUU2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) RETURN
*
      IF( UPPER ) THEN
*
*        Compute the product U * U**H.
*
         DO 10 I = 1, N
            AII = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
               A( I, I ) = AII*AII + DBLE( ZDOTC( N-I, A( I, I+1 ), LDA, A( I, I+1 ), LDA ) )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, DCMPLX( AII ), A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
            ELSE
               CALL ZDSCAL( I, AII, A( 1, I ), 1 )
            END IF
   10    CONTINUE
*
      ELSE
*
*        Compute the product L**H * L.
*
         DO 20 I = 1, N
            AII = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
               A( I, I ) = AII*AII + DBLE( ZDOTC( N-I, A( I+1, I ), 1, A( I+1, I ), 1 ) )
               CALL ZLACGV( I-1, A( I, 1 ), LDA )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, DCMPLX( AII ), A( I, 1 ), LDA )
               CALL ZLACGV( I-1, A( I, 1 ), LDA )
            ELSE
               CALL ZDSCAL( I, AII, A( I, 1 ), LDA )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLAUU2
*
      END
