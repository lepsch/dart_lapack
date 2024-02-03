      SUBROUTINE ZGEQPF( M, N, A, LDA, JPVT, TAU, WORK, RWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             RWORK( * );
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D+0, ONE = 1.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, MA, MN, PVT;
      double             TEMP, TEMP2, TOL3Z;
      COMPLEX*16         AII
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQR2, ZLARF, ZLARFG, ZSWAP, ZUNM2R
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DCMPLX, DCONJG, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DZNRM2;
      // EXTERNAL IDAMAX, DLAMCH, DZNRM2
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGEQPF', -INFO );
         RETURN
      }

      MN = MIN( M, N )
      TOL3Z = SQRT(DLAMCH('Epsilon'))

      // Move initial columns up front

      ITEMP = 1
      for (I = 1; I <= N; I++) { // 10
         if ( JPVT( I ).NE.0 ) {
            if ( I.NE.ITEMP ) {
               zswap(M, A( 1, I ), 1, A( 1, ITEMP ), 1 );
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            } else {
               JPVT( I ) = I
            }
            ITEMP = ITEMP + 1
         } else {
            JPVT( I ) = I
         }
      } // 10
      ITEMP = ITEMP - 1

      // Compute the QR factorization and update remaining columns

      if ( ITEMP.GT.0 ) {
         MA = MIN( ITEMP, M )
         zgeqr2(M, MA, A, LDA, TAU, WORK, INFO );
         if ( MA.LT.N ) {
            zunm2r('Left', 'Conjugate transpose', M, N-MA, MA, A, LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO );
         }
      }

      if ( ITEMP.LT.MN ) {

         // Initialize partial column norms. The first n elements of
         // work store the exact column norms.

         for (I = ITEMP + 1; I <= N; I++) { // 20
            RWORK( I ) = DZNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            RWORK( N+I ) = RWORK( I )
         } // 20

         // Compute factorization

         for (I = ITEMP + 1; I <= MN; I++) { // 40

            // Determine ith pivot column and swap if necessary

            PVT = ( I-1 ) + IDAMAX( N-I+1, RWORK( I ), 1 )

            if ( PVT.NE.I ) {
               zswap(M, A( 1, PVT ), 1, A( 1, I ), 1 );
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               RWORK( PVT ) = RWORK( I )
               RWORK( N+PVT ) = RWORK( N+I )
            }

            // Generate elementary reflector H(i)

            AII = A( I, I )
            zlarfg(M-I+1, AII, A( MIN( I+1, M ), I ), 1, TAU( I ) );
            A( I, I ) = AII

            if ( I.LT.N ) {

               // Apply H(i) to A(i:m,i+1:n) from the left

               AII = A( I, I )
               A( I, I ) = DCMPLX( ONE )
               zlarf('Left', M-I+1, N-I, A( I, I ), 1, DCONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK );
               A( I, I ) = AII
            }

            // Update partial column norms

            for (J = I + 1; J <= N; J++) { // 30
               if ( RWORK( J ).NE.ZERO ) {

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( I, J ) ) / RWORK( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( RWORK( J ) / RWORK( N+J ) )**2
                  if ( TEMP2 .LE. TOL3Z ) {
                     if ( M-I.GT.0 ) {
                        RWORK( J ) = DZNRM2( M-I, A( I+1, J ), 1 )
                        RWORK( N+J ) = RWORK( J )
                     } else {
                        RWORK( J ) = ZERO
                        RWORK( N+J ) = ZERO
                     }
                  } else {
                     RWORK( J ) = RWORK( J )*SQRT( TEMP )
                  }
               }
            } // 30

         } // 40
      }
      RETURN

      // End of ZGEQPF

      }
