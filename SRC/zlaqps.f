      SUBROUTINE ZLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, VN2, AUXV, F, LDF )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KB, LDA, LDF, M, N, NB, OFFSET;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      double             VN1( * ), VN2( * );
      COMPLEX*16         A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      COMPLEX*16         CZERO, CONE
      const              ZERO = 0.0D+0, ONE = 1.0D+0, CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                ITEMP, J, K, LASTRK, LSTICC, PVT, RK;
      double             TEMP, TEMP2, TOL3Z;
      COMPLEX*16         AKK
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM, ZGEMV, ZLARFG, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DCONJG, MAX, MIN, NINT, SQRT
      // ..
      // .. External Functions ..
      int                IDAMAX;
      double             DLAMCH, DZNRM2;
      // EXTERNAL IDAMAX, DLAMCH, DZNRM2
      // ..
      // .. Executable Statements ..

      LASTRK = MIN( M, N+OFFSET )
      LSTICC = 0
      K = 0
      TOL3Z = SQRT(DLAMCH('Epsilon'))

      // Beginning of while loop.

   10 CONTINUE
      if ( ( K.LT.NB ) .AND. ( LSTICC.EQ.0 ) ) {
         K = K + 1
         RK = OFFSET + K

         // Determine ith pivot column and swap if necessary

         PVT = ( K-1 ) + IDAMAX( N-K+1, VN1( K ), 1 )
         if ( PVT.NE.K ) {
            CALL ZSWAP( M, A( 1, PVT ), 1, A( 1, K ), 1 )
            CALL ZSWAP( K-1, F( PVT, 1 ), LDF, F( K, 1 ), LDF )
            ITEMP = JPVT( PVT )
            JPVT( PVT ) = JPVT( K )
            JPVT( K ) = ITEMP
            VN1( PVT ) = VN1( K )
            VN2( PVT ) = VN2( K )
         }

         // Apply previous Householder reflectors to column K:
         // A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**H.

         if ( K.GT.1 ) {
            DO 20 J = 1, K - 1
               F( K, J ) = DCONJG( F( K, J ) )
   20       CONTINUE
            CALL ZGEMV( 'No transpose', M-RK+1, K-1, -CONE, A( RK, 1 ), LDA, F( K, 1 ), LDF, CONE, A( RK, K ), 1 )
            DO 30 J = 1, K - 1
               F( K, J ) = DCONJG( F( K, J ) )
   30       CONTINUE
         }

         // Generate elementary reflector H(k).

         if ( RK.LT.M ) {
            CALL ZLARFG( M-RK+1, A( RK, K ), A( RK+1, K ), 1, TAU( K ) )
         } else {
            CALL ZLARFG( 1, A( RK, K ), A( RK, K ), 1, TAU( K ) )
         }

         AKK = A( RK, K )
         A( RK, K ) = CONE

         // Compute Kth column of F:

         // Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**H*A(RK:M,K).

         if ( K.LT.N ) {
            CALL ZGEMV( 'Conjugate transpose', M-RK+1, N-K, TAU( K ), A( RK, K+1 ), LDA, A( RK, K ), 1, CZERO, F( K+1, K ), 1 )
         }

         // Padding F(1:K,K) with zeros.

         DO 40 J = 1, K
            F( J, K ) = CZERO
   40    CONTINUE

         // Incremental updating of F:
         // F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**H
                     // *A(RK:M,K).

         if ( K.GT.1 ) {
            CALL ZGEMV( 'Conjugate transpose', M-RK+1, K-1, -TAU( K ), A( RK, 1 ), LDA, A( RK, K ), 1, CZERO, AUXV( 1 ), 1 )

            CALL ZGEMV( 'No transpose', N, K-1, CONE, F( 1, 1 ), LDF, AUXV( 1 ), 1, CONE, F( 1, K ), 1 )
         }

         // Update the current row of A:
         // A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**H.

         if ( K.LT.N ) {
            CALL ZGEMM( 'No transpose', 'Conjugate transpose', 1, N-K, K, -CONE, A( RK, 1 ), LDA, F( K+1, 1 ), LDF, CONE, A( RK, K+1 ), LDA )
         }

         // Update partial column norms.

         if ( RK.LT.LASTRK ) {
            DO 50 J = K + 1, N
               if ( VN1( J ).NE.ZERO ) {

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( RK, J ) ) / VN1( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2
                  if ( TEMP2 .LE. TOL3Z ) {
                     VN2( J ) = DBLE( LSTICC )
                     LSTICC = J
                  } else {
                     VN1( J ) = VN1( J )*SQRT( TEMP )
                  }
               }
   50       CONTINUE
         }

         A( RK, K ) = AKK

         // End of while loop.

         GO TO 10
      }
      KB = K
      RK = OFFSET + KB

      // Apply the block reflector to the rest of the matrix:
      // A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
                          // A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**H.

      if ( KB.LT.MIN( N, M-OFFSET ) ) {
         CALL ZGEMM( 'No transpose', 'Conjugate transpose', M-RK, N-KB, KB, -CONE, A( RK+1, 1 ), LDA, F( KB+1, 1 ), LDF, CONE, A( RK+1, KB+1 ), LDA )
      }

      // Recomputation of difficult columns.

   60 CONTINUE
      if ( LSTICC.GT.0 ) {
         ITEMP = NINT( VN2( LSTICC ) )
         VN1( LSTICC ) = DZNRM2( M-RK, A( RK+1, LSTICC ), 1 )

         // NOTE: The computation of VN1( LSTICC ) relies on the fact that
         // SNRM2 does not fail on vectors with norm below the value of
         // SQRT(DLAMCH('S'))

         VN2( LSTICC ) = VN1( LSTICC )
         LSTICC = ITEMP
         GO TO 60
      }

      RETURN

      // End of ZLAQPS

      }
