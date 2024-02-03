      SUBROUTINE SLAQPS( M, N, OFFSET, NB, KB, A, LDA, JPVT, TAU, VN1, VN2, AUXV, F, LDF );

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                KB, LDA, LDF, M, N, NB, OFFSET;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               A( LDA, * ), AUXV( * ), F( LDF, * ), TAU( * ), VN1( * ), VN2( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      int                ITEMP, J, K, LASTRK, LSTICC, PVT, RK;
      REAL               AKK, TEMP, TEMP2, TOL3Z;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMM, SGEMV, SLARFG, SSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, NINT, REAL, SQRT
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH, SNRM2;
      // EXTERNAL ISAMAX, SLAMCH, SNRM2
      // ..
      // .. Executable Statements ..

      LASTRK = min( M, N+OFFSET );
      LSTICC = 0;
      K = 0;
      TOL3Z = sqrt(SLAMCH('Epsilon'));

      // Beginning of while loop.

      } // 10
      if ( ( K < NB ) && ( LSTICC == 0 ) ) {
         K = K + 1;
         RK = OFFSET + K;

         // Determine ith pivot column and swap if necessary

         PVT = ( K-1 ) + ISAMAX( N-K+1, VN1( K ), 1 );
         if ( PVT != K ) {
            sswap(M, A( 1, PVT ), 1, A( 1, K ), 1 );
            sswap(K-1, F( PVT, 1 ), LDF, F( K, 1 ), LDF );
            ITEMP = JPVT( PVT );
            JPVT( PVT ) = JPVT( K );
            JPVT( K ) = ITEMP;
            VN1( PVT ) = VN1( K );
            VN2( PVT ) = VN2( K );
         }

         // Apply previous Householder reflectors to column K:
         // A(RK:M,K) := A(RK:M,K) - A(RK:M,1:K-1)*F(K,1:K-1)**T.

         if ( K > 1 ) {
            sgemv('No transpose', M-RK+1, K-1, -ONE, A( RK, 1 ), LDA, F( K, 1 ), LDF, ONE, A( RK, K ), 1 );
         }

         // Generate elementary reflector H(k).

         if ( RK < M ) {
            slarfg(M-RK+1, A( RK, K ), A( RK+1, K ), 1, TAU( K ) );
         } else {
            slarfg(1, A( RK, K ), A( RK, K ), 1, TAU( K ) );
         }

         AKK = A( RK, K );
         A( RK, K ) = ONE;

         // Compute Kth column of F:

         // Compute  F(K+1:N,K) := tau(K)*A(RK:M,K+1:N)**T*A(RK:M,K).

         if ( K < N ) {
            sgemv('Transpose', M-RK+1, N-K, TAU( K ), A( RK, K+1 ), LDA, A( RK, K ), 1, ZERO, F( K+1, K ), 1 );
         }

         // Padding F(1:K,K) with zeros.

         for (J = 1; J <= K; J++) { // 20
            F( J, K ) = ZERO;
         } // 20

         // Incremental updating of F:
         // F(1:N,K) := F(1:N,K) - tau(K)*F(1:N,1:K-1)*A(RK:M,1:K-1)**T
                     // *A(RK:M,K).

         if ( K > 1 ) {
            sgemv('Transpose', M-RK+1, K-1, -TAU( K ), A( RK, 1 ), LDA, A( RK, K ), 1, ZERO, AUXV( 1 ), 1 );

            sgemv('No transpose', N, K-1, ONE, F( 1, 1 ), LDF, AUXV( 1 ), 1, ONE, F( 1, K ), 1 );
         }

         // Update the current row of A:
         // A(RK,K+1:N) := A(RK,K+1:N) - A(RK,1:K)*F(K+1:N,1:K)**T.

         if ( K < N ) {
            sgemv('No transpose', N-K, K, -ONE, F( K+1, 1 ), LDF, A( RK, 1 ), LDA, ONE, A( RK, K+1 ), LDA );
         }

         // Update partial column norms.

         if ( RK < LASTRK ) {
            for (J = K + 1; J <= N; J++) { // 30
               if ( VN1( J ) != ZERO ) {

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( RK, J ) ) / VN1( J );
                  TEMP = max( ZERO, ( ONE+TEMP )*( ONE-TEMP ) );
                  TEMP2 = TEMP*( VN1( J ) / VN2( J ) )**2;
                  if ( TEMP2 <= TOL3Z ) {
                     VN2( J ) = REAL( LSTICC );
                     LSTICC = J;
                  } else {
                     VN1( J ) = VN1( J )*sqrt( TEMP );
                  }
               }
            } // 30
         }

         A( RK, K ) = AKK;

         // End of while loop.

         GO TO 10;
      }
      KB = K;
      RK = OFFSET + KB;

      // Apply the block reflector to the rest of the matrix:
      // A(OFFSET+KB+1:M,KB+1:N) := A(OFFSET+KB+1:M,KB+1:N) -
                          // A(OFFSET+KB+1:M,1:KB)*F(KB+1:N,1:KB)**T.

      if ( KB < min( N, M-OFFSET ) ) {
         sgemm('No transpose', 'Transpose', M-RK, N-KB, KB, -ONE, A( RK+1, 1 ), LDA, F( KB+1, 1 ), LDF, ONE, A( RK+1, KB+1 ), LDA );
      }

      // Recomputation of difficult columns.

      } // 40
      if ( LSTICC > 0 ) {
         ITEMP = NINT( VN2( LSTICC ) );
         VN1( LSTICC ) = SNRM2( M-RK, A( RK+1, LSTICC ), 1 );

         // NOTE: The computation of VN1( LSTICC ) relies on the fact that
         // SNRM2 does not fail on vectors with norm below the value of
         // sqrt(DLAMCH('S'))

         VN2( LSTICC ) = VN1( LSTICC );
         LSTICC = ITEMP;
         GO TO 40;
      }

      return;

      // End of SLAQPS

      }
