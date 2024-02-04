      void cggsvp3(JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, RWORK, TAU, WORK, LWORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // IMPLICIT NONE

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P, LWORK;
      double               TOLA, TOLB;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      double               RWORK( * );
      Complex            A( LDA, * ), B( LDB, * ), Q( LDQ, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               FORWRD, WANTQ, WANTU, WANTV, LQUERY;
      int                I, J, LWKOPT;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEQP3, CGEQR2, CGERQ2, CLACPY, CLAPMT, CLASET, CUNG2R, CUNM2R, CUNMR2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      WANTU = lsame( JOBU, 'U' );
      WANTV = lsame( JOBV, 'V' );
      WANTQ = lsame( JOBQ, 'Q' );
      FORWRD = true;
      LQUERY = ( LWORK == -1 );
      LWKOPT = 1;

      // Test the input arguments

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
      } else if ( LWORK < 1 && !LQUERY ) {
         INFO = -24;
      }

      // Compute workspace

      if ( INFO == 0 ) {
         cgeqp3(P, N, B, LDB, IWORK, TAU, WORK, -1, RWORK, INFO );
         LWKOPT = INT( WORK ( 1 ) );
         if ( WANTV ) {
            LWKOPT = max( LWKOPT, P );
         }
         LWKOPT = max( LWKOPT, min( N, P ) );
         LWKOPT = max( LWKOPT, M );
         if ( WANTQ ) {
            LWKOPT = max( LWKOPT, N );
         }
         cgeqp3(M, N, A, LDA, IWORK, TAU, WORK, -1, RWORK, INFO );
         LWKOPT = max( LWKOPT, INT( WORK ( 1 ) ) );
         LWKOPT = max( 1, LWKOPT );
         WORK[1] = CMPLX( LWKOPT );
      }

      if ( INFO != 0 ) {
         xerbla('CGGSVP3', -INFO );
         return;
      }
      if ( LQUERY ) {
         return;
      }

      // QR with column pivoting of B: B*P = V*( S11 S12 )
                                            // (  0   0  )

      for (I = 1; I <= N; I++) { // 10
         IWORK[I] = 0;
      } // 10
      cgeqp3(P, N, B, LDB, IWORK, TAU, WORK, LWORK, RWORK, INFO );

      // Update A := A*P

      clapmt(FORWRD, M, N, A, LDA, IWORK );

      // Determine the effective rank of matrix B.

      L = 0;
      for (I = 1; I <= min( P, N ); I++) { // 20
         if( ( B( I, I ) ).abs() > TOLB ) L = L + 1;
      } // 20

      if ( WANTV ) {

         // Copy the details of V, and form V.

         claset('Full', P, P, CZERO, CZERO, V, LDV );
         if (P > 1) clacpy( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ), LDV );
         cung2r(P, P, min( P, N ), V, LDV, TAU, WORK, INFO );
      }

      // Clean up B

      for (J = 1; J <= L - 1; J++) { // 40
         for (I = J + 1; I <= L; I++) { // 30
            B[I, J] = CZERO;
         } // 30
      } // 40
      if (P > L) claset( 'Full', P-L, N, CZERO, CZERO, B( L+1, 1 ), LDB );

      if ( WANTQ ) {

         // Set Q = I and Update Q := Q*P

         claset('Full', N, N, CZERO, CONE, Q, LDQ );
         clapmt(FORWRD, N, N, Q, LDQ, IWORK );
      }

      if ( P >= L && N != L ) {

         // RQ factorization of ( S11 S12 ) = ( 0 S12 )*Z

         cgerq2(L, N, B, LDB, TAU, WORK, INFO );

         // Update A := A*Z**H

         cunmr2('Right', 'Conjugate transpose', M, N, L, B, LDB, TAU, A, LDA, WORK, INFO );
         if ( WANTQ ) {

            // Update Q := Q*Z**H

            cunmr2('Right', 'Conjugate transpose', N, N, L, B, LDB, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up B

         claset('Full', L, N-L, CZERO, CZERO, B, LDB );
         for (J = N - L + 1; J <= N; J++) { // 60
            for (I = J - N + L + 1; I <= L; I++) { // 50
               B[I, J] = CZERO;
            } // 50
         } // 60

      }

      // Let              N-L     L
                 // A = ( A11    A12 ) M,

      // then the following does the complete QR decomposition of A11:

               // A11 = U*(  0  T12 )*P1**H
                       // (  0   0  )

      for (I = 1; I <= N - L; I++) { // 70
         IWORK[I] = 0;
      } // 70
      cgeqp3(M, N-L, A, LDA, IWORK, TAU, WORK, LWORK, RWORK, INFO );

      // Determine the effective rank of A11

      K = 0;
      for (I = 1; I <= min( M, N-L ); I++) { // 80
         if( ( A( I, I ) ).abs() > TOLA ) K = K + 1;
      } // 80

      // Update A12 := U**H*A12, where A12 = A( 1:M, N-L+1:N )

      cunm2r('Left', 'Conjugate transpose', M, L, min( M, N-L ), A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO );

      if ( WANTU ) {

         // Copy the details of U, and form U

         claset('Full', M, M, CZERO, CZERO, U, LDU );
         if (M > 1) clacpy( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ), LDU );
         cung2r(M, M, min( M, N-L ), U, LDU, TAU, WORK, INFO );
      }

      if ( WANTQ ) {

         // Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1

         clapmt(FORWRD, N, N-L, Q, LDQ, IWORK );
      }

      // Clean up A: set the strictly lower triangular part of
      // A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.

      for (J = 1; J <= K - 1; J++) { // 100
         for (I = J + 1; I <= K; I++) { // 90
            A[I, J] = CZERO;
         } // 90
      } // 100
      if (M > K) claset( 'Full', M-K, N-L, CZERO, CZERO, A( K+1, 1 ), LDA );

      if ( N-L > K ) {

         // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

         cgerq2(K, N-L, A, LDA, TAU, WORK, INFO );

         if ( WANTQ ) {

            // Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**H

            cunmr2('Right', 'Conjugate transpose', N, N-L, K, A, LDA, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up A

         claset('Full', K, N-L-K, CZERO, CZERO, A, LDA );
         for (J = N - L - K + 1; J <= N - L; J++) { // 120
            for (I = J - N + L + K + 1; I <= K; I++) { // 110
               A[I, J] = CZERO;
            } // 110
         } // 120

      }

      if ( M > K ) {

         // QR factorization of A( K+1:M,N-L+1:N )

         cgeqr2(M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO );

         if ( WANTU ) {

            // Update U(:,K+1:M) := U(:,K+1:M)*U1

            cunm2r('Right', 'No transpose', M, M-K, min( M-K, L ), A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU, WORK, INFO );
         }

         // Clean up

         for (J = N - L + 1; J <= N; J++) { // 140
            for (I = J - N + K + L + 1; I <= M; I++) { // 130
               A[I, J] = CZERO;
            } // 130
         } // 140

      }

      WORK[1] = CMPLX( LWKOPT );
      return;
      }