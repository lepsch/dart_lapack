      void zggsvp(final int JOBU, final int JOBV, final int JOBQ, final int M, final int P, final int N, final Matrix<double> A_, final int LDA, final Matrix<double> B_, final int LDB, final int TOLA, final int TOLB, final int K, final int L, final Matrix<double> U_, final int LDU, final Matrix<double> V_, final int LDV, final Matrix<double> Q_, final int LDQ, final Array<int> IWORK_, final Array<double> RWORK_, final int TAU, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.having();
  final B = B_.having();
  final U = U_.having();
  final V = V_.having();
  final Q = Q_.having();
  final IWORK = IWORK_.having();
  final RWORK = RWORK_.having();
  final _WORK = _WORK_.having();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P;
      double             TOLA, TOLB;
      int                IWORK( * );
      double             RWORK( * );
      Complex         A( LDA, * ), B( LDB, * ), Q( LDQ, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

      Complex         CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               FORWRD, WANTQ, WANTU, WANTV;
      int                I, J;
      Complex         T;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQPF, ZGEQR2, ZGERQ2, ZLACPY, ZLAPMT, ZLASET, ZUNG2R, ZUNM2R, ZUNMR2
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, DBLE, DIMAG, MAX, MIN
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1[T] = ( T.toDouble() ).abs() + ( DIMAG( T ) ).abs();

      // Test the input parameters

      WANTU = lsame( JOBU, 'U' );
      WANTV = lsame( JOBV, 'V' );
      WANTQ = lsame( JOBQ, 'Q' );
      FORWRD = true;

      INFO = 0;
      if ( !( WANTU || lsame( JOBU, 'N' ) ) ) {
         INFO = -1;
      } else if ( !( WANTV || lsame( JOBV, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( WANTQ || lsame( JOBQ, 'N' ) ) ) {
         INFO = -3;
      } else if ( M < 0 ) {
         INFO = -4;
      } else if ( P < 0 ) {
         INFO = -5;
      } else if ( N < 0 ) {
         INFO = -6;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -8;
      } else if ( LDB < max( 1, P ) ) {
         INFO = -10;
      } else if ( LDU < 1 || ( WANTU && LDU < M ) ) {
         INFO = -16;
      } else if ( LDV < 1 || ( WANTV && LDV < P ) ) {
         INFO = -18;
      } else if ( LDQ < 1 || ( WANTQ && LDQ < N ) ) {
         INFO = -20;
      }
      if ( INFO != 0 ) {
         xerbla('ZGGSVP', -INFO );
         return;
      }

      // QR with column pivoting of B: B*P = V*( S11 S12 )
      //                                       (  0   0  )

      for (I = 1; I <= N; I++) { // 10
         IWORK[I] = 0;
      } // 10
      zgeqpf(P, N, B, LDB, IWORK, TAU, WORK, RWORK, INFO );

      // Update A := A*P

      zlapmt(FORWRD, M, N, A, LDA, IWORK );

      // Determine the effective rank of matrix B.

      L = 0;
      for (I = 1; I <= min( P, N ); I++) { // 20
         if( CABS1( B( I, I ) ) > TOLB ) L++;
      } // 20

      if ( WANTV ) {

         // Copy the details of V, and form V.

         zlaset('Full', P, P, CZERO, CZERO, V, LDV );
         if (P > 1) zlacpy( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ), LDV );
         zung2r(P, P, min( P, N ), V, LDV, TAU, WORK, INFO );
      }

      // Clean up B

      for (J = 1; J <= L - 1; J++) { // 40
         for (I = J + 1; I <= L; I++) { // 30
            B[I][J] = CZERO;
         } // 30
      } // 40
      if (P > L) zlaset( 'Full', P-L, N, CZERO, CZERO, B( L+1, 1 ), LDB );

      if ( WANTQ ) {

         // Set Q = I and Update Q := Q*P

         zlaset('Full', N, N, CZERO, CONE, Q, LDQ );
         zlapmt(FORWRD, N, N, Q, LDQ, IWORK );
      }

      if ( P >= L && N != L ) {

         // RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z

         zgerq2(L, N, B, LDB, TAU, WORK, INFO );

         // Update A := A*Z**H

         zunmr2('Right', 'Conjugate transpose', M, N, L, B, LDB, TAU, A, LDA, WORK, INFO );
         if ( WANTQ ) {

            // Update Q := Q*Z**H

            zunmr2('Right', 'Conjugate transpose', N, N, L, B, LDB, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up B

         zlaset('Full', L, N-L, CZERO, CZERO, B, LDB );
         for (J = N - L + 1; J <= N; J++) { // 60
            for (I = J - N + L + 1; I <= L; I++) { // 50
               B[I][J] = CZERO;
            } // 50
         } // 60

      }

      // Let              N-L     L
      //            A = ( A11    A12 ) M,

      // then the following does the complete QR decomposition of A11:

               // A11 = U*(  0  T12 )*P1**H
               //         (  0   0  )

      for (I = 1; I <= N - L; I++) { // 70
         IWORK[I] = 0;
      } // 70
      zgeqpf(M, N-L, A, LDA, IWORK, TAU, WORK, RWORK, INFO );

      // Determine the effective rank of A11

      K = 0;
      for (I = 1; I <= min( M, N-L ); I++) { // 80
         if( CABS1( A( I, I ) ) > TOLA ) K++;
      } // 80

      // Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N )

      zunm2r('Left', 'Conjugate transpose', M, L, min( M, N-L ), A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO );

      if ( WANTU ) {

         // Copy the details of U, and form U

         zlaset('Full', M, M, CZERO, CZERO, U, LDU );
         if (M > 1) zlacpy( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ), LDU );
         zung2r(M, M, min( M, N-L ), U, LDU, TAU, WORK, INFO );
      }

      if ( WANTQ ) {

         // Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1

         zlapmt(FORWRD, N, N-L, Q, LDQ, IWORK );
      }

      // Clean up A: set the strictly lower triangular part of
      // A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.

      for (J = 1; J <= K - 1; J++) { // 100
         for (I = J + 1; I <= K; I++) { // 90
            A[I][J] = CZERO;
         } // 90
      } // 100
      if (M > K) zlaset( 'Full', M-K, N-L, CZERO, CZERO, A( K+1, 1 ), LDA );

      if ( N-L > K ) {

         // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

         zgerq2(K, N-L, A, LDA, TAU, WORK, INFO );

         if ( WANTQ ) {

            // Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H

            zunmr2('Right', 'Conjugate transpose', N, N-L, K, A, LDA, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up A

         zlaset('Full', K, N-L-K, CZERO, CZERO, A, LDA );
         for (J = N - L - K + 1; J <= N - L; J++) { // 120
            for (I = J - N + L + K + 1; I <= K; I++) { // 110
               A[I][J] = CZERO;
            } // 110
         } // 120

      }

      if ( M > K ) {

         // QR factorization of A( K+1:M,N-L+1:N )

         zgeqr2(M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO );

         if ( WANTU ) {

            // Update U(:,K+1:M) := U(:,K+1:M)*U1

            zunm2r('Right', 'No transpose', M, M-K, min( M-K, L ), A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU, WORK, INFO );
         }

         // Clean up

         for (J = N - L + 1; J <= N; J++) { // 140
            for (I = J - N + K + L + 1; I <= M; I++) { // 130
               A[I][J] = CZERO;
            } // 130
         } // 140

      }

      }
