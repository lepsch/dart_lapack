      SUBROUTINE SGGSVP( JOBU, JOBV, JOBQ, M, P, N, A, LDA, B, LDB, TOLA, TOLB, K, L, U, LDU, V, LDV, Q, LDQ, IWORK, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBQ, JOBU, JOBV;
      int                INFO, K, L, LDA, LDB, LDQ, LDU, LDV, M, N, P;
      REAL               TOLA, TOLB
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), Q( LDQ, * ), TAU( * ), U( LDU, * ), V( LDV, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               FORWRD, WANTQ, WANTU, WANTV;
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQPF, SGEQR2, SGERQ2, SLACPY, SLAPMT, SLASET, SORG2R, SORM2R, SORMR2, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      WANTU = LSAME( JOBU, 'U' )
      WANTV = LSAME( JOBV, 'V' )
      WANTQ = LSAME( JOBQ, 'Q' )
      FORWRD = .TRUE.

      INFO = 0
      if ( .NOT.( WANTU .OR. LSAME( JOBU, 'N' ) ) ) {
         INFO = -1
      } else if ( .NOT.( WANTV .OR. LSAME( JOBV, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( WANTQ .OR. LSAME( JOBQ, 'N' ) ) ) {
         INFO = -3
      } else if ( M.LT.0 ) {
         INFO = -4
      } else if ( P.LT.0 ) {
         INFO = -5
      } else if ( N.LT.0 ) {
         INFO = -6
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -8
      } else if ( LDB.LT.MAX( 1, P ) ) {
         INFO = -10
      } else if ( LDU.LT.1 .OR. ( WANTU .AND. LDU.LT.M ) ) {
         INFO = -16
      } else if ( LDV.LT.1 .OR. ( WANTV .AND. LDV.LT.P ) ) {
         INFO = -18
      } else if ( LDQ.LT.1 .OR. ( WANTQ .AND. LDQ.LT.N ) ) {
         INFO = -20
      }
      if ( INFO.NE.0 ) {
         xerbla('SGGSVP', -INFO );
         RETURN
      }

      // QR with column pivoting of B: B*P = V*( S11 S12 )
                                            // (  0   0  )

      for (I = 1; I <= N; I++) { // 10
         IWORK( I ) = 0
      } // 10
      sgeqpf(P, N, B, LDB, IWORK, TAU, WORK, INFO );

      // Update A := A*P

      slapmt(FORWRD, M, N, A, LDA, IWORK );

      // Determine the effective rank of matrix B.

      L = 0
      DO 20 I = 1, MIN( P, N )
         IF( ABS( B( I, I ) ).GT.TOLB ) L = L + 1
      } // 20

      if ( WANTV ) {

         // Copy the details of V, and form V.

         slaset('Full', P, P, ZERO, ZERO, V, LDV );
         IF( P.GT.1 ) CALL SLACPY( 'Lower', P-1, N, B( 2, 1 ), LDB, V( 2, 1 ), LDV )
         sorg2r(P, P, MIN( P, N ), V, LDV, TAU, WORK, INFO );
      }

      // Clean up B

      for (J = 1; J <= L - 1; J++) { // 40
         for (I = J + 1; I <= L; I++) { // 30
            B( I, J ) = ZERO
         } // 30
      } // 40
      IF( P.GT.L ) CALL SLASET( 'Full', P-L, N, ZERO, ZERO, B( L+1, 1 ), LDB )

      if ( WANTQ ) {

         // Set Q = I and Update Q := Q*P

         slaset('Full', N, N, ZERO, ONE, Q, LDQ );
         slapmt(FORWRD, N, N, Q, LDQ, IWORK );
      }

      if ( P.GE.L .AND. N.NE.L ) {

         // RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z

         sgerq2(L, N, B, LDB, TAU, WORK, INFO );

         // Update A := A*Z**T

         sormr2('Right', 'Transpose', M, N, L, B, LDB, TAU, A, LDA, WORK, INFO );

         if ( WANTQ ) {

            // Update Q := Q*Z**T

            sormr2('Right', 'Transpose', N, N, L, B, LDB, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up B

         slaset('Full', L, N-L, ZERO, ZERO, B, LDB );
         for (J = N - L + 1; J <= N; J++) { // 60
            for (I = J - N + L + 1; I <= L; I++) { // 50
               B( I, J ) = ZERO
            } // 50
         } // 60

      }

      // Let              N-L     L
                 // A = ( A11    A12 ) M,

      // then the following does the complete QR decomposition of A11:

               // A11 = U*(  0  T12 )*P1**T
                       // (  0   0  )

      for (I = 1; I <= N - L; I++) { // 70
         IWORK( I ) = 0
      } // 70
      sgeqpf(M, N-L, A, LDA, IWORK, TAU, WORK, INFO );

      // Determine the effective rank of A11

      K = 0
      DO 80 I = 1, MIN( M, N-L )
         IF( ABS( A( I, I ) ).GT.TOLA ) K = K + 1
      } // 80

      // Update A12 := U**T*A12, where A12 = A( 1:M, N-L+1:N )

      sorm2r('Left', 'Transpose', M, L, MIN( M, N-L ), A, LDA, TAU, A( 1, N-L+1 ), LDA, WORK, INFO );

      if ( WANTU ) {

         // Copy the details of U, and form U

         slaset('Full', M, M, ZERO, ZERO, U, LDU );
         IF( M.GT.1 ) CALL SLACPY( 'Lower', M-1, N-L, A( 2, 1 ), LDA, U( 2, 1 ), LDU )
         sorg2r(M, M, MIN( M, N-L ), U, LDU, TAU, WORK, INFO );
      }

      if ( WANTQ ) {

         // Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1

         slapmt(FORWRD, N, N-L, Q, LDQ, IWORK );
      }

      // Clean up A: set the strictly lower triangular part of
      // A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0.

      for (J = 1; J <= K - 1; J++) { // 100
         for (I = J + 1; I <= K; I++) { // 90
            A( I, J ) = ZERO
         } // 90
      } // 100
      IF( M.GT.K ) CALL SLASET( 'Full', M-K, N-L, ZERO, ZERO, A( K+1, 1 ), LDA )

      if ( N-L.GT.K ) {

         // RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1

         sgerq2(K, N-L, A, LDA, TAU, WORK, INFO );

         if ( WANTQ ) {

            // Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1**T

            sormr2('Right', 'Transpose', N, N-L, K, A, LDA, TAU, Q, LDQ, WORK, INFO );
         }

         // Clean up A

         slaset('Full', K, N-L-K, ZERO, ZERO, A, LDA );
         for (J = N - L - K + 1; J <= N - L; J++) { // 120
            for (I = J - N + L + K + 1; I <= K; I++) { // 110
               A( I, J ) = ZERO
            } // 110
         } // 120

      }

      if ( M.GT.K ) {

         // QR factorization of A( K+1:M,N-L+1:N )

         sgeqr2(M-K, L, A( K+1, N-L+1 ), LDA, TAU, WORK, INFO );

         if ( WANTU ) {

            // Update U(:,K+1:M) := U(:,K+1:M)*U1

            sorm2r('Right', 'No transpose', M, M-K, MIN( M-K, L ), A( K+1, N-L+1 ), LDA, TAU, U( 1, K+1 ), LDU, WORK, INFO );
         }

         // Clean up

         for (J = N - L + 1; J <= N; J++) { // 140
            for (I = J - N + K + L + 1; I <= M; I++) { // 130
               A( I, J ) = ZERO
            } // 130
         } // 140

      }

      RETURN

      // End of SGGSVP

      }
