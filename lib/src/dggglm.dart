      void dggglm(N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), X( * ), Y( * );
      // ..

// ===================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0, ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DGEMV, DGGQRF, DORMQR, DORMRQ, DTRTRS, XERBLA
      // ..
      // .. External Functions ..
      //- int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

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
            NB1 = ILAENV( 1, 'DGEQRF', ' ', N, M, -1, -1 );
            NB2 = ILAENV( 1, 'DGERQF', ' ', N, M, -1, -1 );
            NB3 = ILAENV( 1, 'DORMQR', ' ', N, M, P, -1 );
            NB4 = ILAENV( 1, 'DORMRQ', ' ', N, M, P, -1 );
            NB = max( NB1, NB2, NB3, NB4 );
            LWKMIN = M + N + P;
            LWKOPT = M + NP + max( N, P )*NB;
         }
         WORK[1] = LWKOPT;

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGGGLM', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if ( N == 0 ) {
         for (I = 1; I <= M; I++) {
            X[I] = ZERO;
         }
         for (I = 1; I <= P; I++) {
            Y[I] = ZERO;
         }
         return;
      }

      // Compute the GQR factorization of matrices A and B:

           // Q**T*A = ( R11 ) M,    Q**T*B*Z**T = ( T11   T12 ) M
                    // (  0  ) N-M                 (  0    T22 ) N-M
                       // M                         M+P-N  N-M

      // where R11 and T22 are upper triangular, and Q and Z are
      // orthogonal.

      dggqrf(N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = INT( WORK( M+NP+1 ) );

      // Update left-hand-side vector d = Q**T*d = ( d1 ) M
                                                // ( d2 ) N-M

      dormqr('Left', 'Transpose', N, 1, M, A, LDA, WORK, D, max( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = max( LOPT, INT( WORK( M+NP+1 ) ) );

      // Solve T22*y2 = d2 for y2

      if ( N > M ) {
         dtrtrs('Upper', 'No transpose', 'Non unit', N-M, 1, B( M+1, M+P-N+1 ), LDB, D( M+1 ), N-M, INFO );

         if ( INFO > 0 ) {
            INFO = 1;
            return;
         }

         dcopy(N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 );
      }

      // Set y1 = 0

      for (I = 1; I <= M + P - N; I++) { // 10
         Y[I] = ZERO;
      } // 10

      // Update d1 = d1 - T12*y2

      dgemv('No transpose', M, N-M, -ONE, B( 1, M+P-N+1 ), LDB, Y( M+P-N+1 ), 1, ONE, D, 1 );

      // Solve triangular system: R11*x = d1

      if ( M > 0 ) {
         dtrtrs('Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D, M, INFO );

         if ( INFO > 0 ) {
            INFO = 2;
            return;
         }

         // Copy D to X

         dcopy(M, D, 1, X, 1 );
      }

      // Backward transformation y = Z**T *y

      dormrq('Left', 'Transpose', P, 1, NP, B( max( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y, max( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      WORK[1] = M + NP + max( LOPT, INT( WORK( M+NP+1 ) ) );

      return;
      }