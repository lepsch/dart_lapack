      SUBROUTINE ZLARGE( N, A, LDA, ISEED, WORK, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             WN;
      COMPLEX*16         TAU, WA, WB
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMV, ZGERC, ZLARNV, ZSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, MAX
      // ..
      // .. External Functions ..
      double             DZNRM2;
      // EXTERNAL DZNRM2
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -3
      }
      if ( INFO.LT.0 ) {
         CALL XERBLA( 'ZLARGE', -INFO )
         RETURN
      }

      // pre- and post-multiply A by random unitary matrix

      DO 10 I = N, 1, -1

         // generate random reflection

         CALL ZLARNV( 3, ISEED, N-I+1, WORK )
         WN = DZNRM2( N-I+1, WORK, 1 )
         WA = ( WN / ABS( WORK( 1 ) ) )*WORK( 1 )
         if ( WN.EQ.ZERO ) {
            TAU = ZERO
         } else {
            WB = WORK( 1 ) + WA
            CALL ZSCAL( N-I, ONE / WB, WORK( 2 ), 1 )
            WORK( 1 ) = ONE
            TAU = DBLE( WB / WA )
         }

         // multiply A(i:n,1:n) by random reflection from the left

         CALL ZGEMV( 'Conjugate transpose', N-I+1, N, ONE, A( I, 1 ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )          CALL ZGERC( N-I+1, N, -TAU, WORK, 1, WORK( N+1 ), 1, A( I, 1 ), LDA )

         // multiply A(1:n,i:n) by random reflection from the right

         CALL ZGEMV( 'No transpose', N, N-I+1, ONE, A( 1, I ), LDA, WORK, 1, ZERO, WORK( N+1 ), 1 )          CALL ZGERC( N, N-I+1, -TAU, WORK( N+1 ), 1, WORK, 1, A( 1, I ), LDA )
   10 CONTINUE
      RETURN

      // End of ZLARGE

      }
