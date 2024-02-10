      void cggglm(N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, Box<int> INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, LDA, LDB, LWORK, M, N, P;
      Complex            A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), X( * ), Y( * );
      // ..

// ===================================================================

      // .. Parameters ..
      Complex            CZERO, CONE;
      const              CZERO = ( 0.0, 0.0 ), CONE = ( 1.0, 0.0 ) ;
      bool               LQUERY;
      int                I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMV, CGGQRF, CTRTRS, CUNMQR, CUNMRQ, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      //- REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN

      // Test the input parameters

      INFO = 0;
      NP = min( N, P );
      LQUERY = ( LWORK == -1 );
      if ( N < 0 ) {
         INFO = -1;
      } else if ( M < 0 || M > N ) {
         INFO = -2;
      } else if ( P < 0 || P < N-M ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, N ) ) {
         INFO = -7;
      }

      // Calculate workspace

      if ( INFO == 0) {
         if ( N == 0 ) {
            LWKMIN = 1;
            LWKOPT = 1;
         } else {
            NB1 = ilaenv( 1, 'CGEQRF', ' ', N, M, -1, -1 );
            NB2 = ilaenv( 1, 'CGERQF', ' ', N, M, -1, -1 );
            NB3 = ilaenv( 1, 'CUNMQR', ' ', N, M, P, -1 );
            NB4 = ilaenv( 1, 'CUNMRQ', ' ', N, M, P, -1 );
            NB = max( NB1, NB2, NB3, NB4 );
            LWKMIN = M + N + P;
            LWKOPT = M + NP + max( N, P )*NB;
         }
         WORK[1] = SROUNDUP_LWORK(LWKOPT);

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGGGLM', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         for (I = 1; I <= M; I++) {
            X[I] = CZERO;
         }
         for (I = 1; I <= P; I++) {
            Y[I] = CZERO;
         }
         return;
      }

      // Compute the GQR factorization of matrices A and B:

           // Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M
           //          (  0  ) N-M                 (  0    T22 ) N-M
           //             M                         M+P-N  N-M

      // where R11 and T22 are upper triangular, and Q and Z are
      // unitary.

      cggqrf(N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = INT( WORK( M+NP+1 ) );

      // Update left-hand-side vector d = Q**H*d = ( d1 ) M
      //                                           ( d2 ) N-M

      cunmqr('Left', 'Conjugate transpose', N, 1, M, A, LDA, WORK, D, max( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = max( LOPT, INT( WORK( M+NP+1 ) ) );

      // Solve T22*y2 = d2 for y2

      if ( N > M ) {
         ctrtrs('Upper', 'No transpose', 'Non unit', N-M, 1, B( M+1, M+P-N+1 ), LDB, D( M+1 ), N-M, INFO );

         if ( INFO > 0 ) {
            INFO = 1;
            return;
         }

         ccopy(N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 );
      }

      // Set y1 = 0

      for (I = 1; I <= M + P - N; I++) { // 10
         Y[I] = CZERO;
      } // 10

      // Update d1 = d1 - T12*y2

      cgemv('No transpose', M, N-M, -CONE, B( 1, M+P-N+1 ), LDB, Y( M+P-N+1 ), 1, CONE, D, 1 );

      // Solve triangular system: R11*x = d1

      if ( M > 0 ) {
         ctrtrs('Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D, M, INFO );

         if ( INFO > 0 ) {
            INFO = 2;
            return;
         }

         // Copy D to X

         ccopy(M, D, 1, X, 1 );
      }

      // Backward transformation y = Z**H *y

      cunmrq('Left', 'Conjugate transpose', P, 1, NP, B( max( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y, max( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      WORK[1] = M + NP + max( LOPT, INT( WORK( M+NP+1 ) ) );

      }
