      SUBROUTINE ZLAUU2( UPLO, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE
      const              ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               UPPER;
      int                I;
      double             AII;
      // ..
      // .. External Functions ..
      bool               LSAME;
      COMPLEX*16         ZDOTC
      // EXTERNAL LSAME, ZDOTC
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZDSCAL, ZGEMV, ZLACGV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCMPLX, MAX
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'ZLAUU2', -INFO )
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      if ( UPPER ) {

         // Compute the product U * U**H.

         DO 10 I = 1, N
            AII = DBLE( A( I, I ) )
            if ( I.LT.N ) {
               A( I, I ) = AII*AII + DBLE( ZDOTC( N-I, A( I, I+1 ), LDA, A( I, I+1 ), LDA ) )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ), LDA, A( I, I+1 ), LDA, DCMPLX( AII ), A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
            } else {
               CALL ZDSCAL( I, AII, A( 1, I ), 1 )
            }
   10    CONTINUE

      } else {

         // Compute the product L**H * L.

         DO 20 I = 1, N
            AII = DBLE( A( I, I ) )
            if ( I.LT.N ) {
               A( I, I ) = AII*AII + DBLE( ZDOTC( N-I, A( I+1, I ), 1, A( I+1, I ), 1 ) )
               CALL ZLACGV( I-1, A( I, 1 ), LDA )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE, A( I+1, 1 ), LDA, A( I+1, I ), 1, DCMPLX( AII ), A( I, 1 ), LDA )
               CALL ZLACGV( I-1, A( I, 1 ), LDA )
            } else {
               CALL ZDSCAL( I, AII, A( I, 1 ), LDA )
            }
   20    CONTINUE
      }

      RETURN

      // End of ZLAUU2

      }
