      void dgglse(M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), B( LDB, * ), C( * ), D( * ), WORK( * ), X( * );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE;
      const              ONE = 1.0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4, NR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DGGRQF, DORMQR, DORMRQ, DTRMV, DTRTRS, XERBLA
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
      MN = min( M, N );
      LQUERY = ( LWORK == -1 );
      if ( M < 0 ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( P < 0 || P > N || P < N-M ) {
         INFO = -3;
      } else if ( LDA < max( 1, M ) ) {
         INFO = -5;
      } else if ( LDB < max( 1, P ) ) {
         INFO = -7;
      }

      // Calculate workspace

      if ( INFO == 0) {
         if ( N == 0 ) {
            LWKMIN = 1;
            LWKOPT = 1;
         } else {
            NB1 = ILAENV( 1, 'DGEQRF', ' ', M, N, -1, -1 );
            NB2 = ILAENV( 1, 'DGERQF', ' ', M, N, -1, -1 );
            NB3 = ILAENV( 1, 'DORMQR', ' ', M, N, P, -1 );
            NB4 = ILAENV( 1, 'DORMRQ', ' ', M, N, P, -1 );
            NB = max( NB1, NB2, NB3, NB4 );
            LWKMIN = M + N + P;
            LWKOPT = P + MN + max( M, N )*NB;
         }
         WORK[1] = LWKOPT;

         if ( LWORK < LWKMIN && !LQUERY ) {
            INFO = -12;
         }
      }

      if ( INFO != 0 ) {
         xerbla('DGGLSE', -INFO );
         return;
      } else if ( LQUERY ) {
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Compute the GRQ factorization of matrices B and A:

             // B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P
                         // N-P  P                     (  0  R22 ) M+P-N
                                                       // N-P  P

      // where T12 and R11 are upper triangular, and Q and Z are
      // orthogonal.

      dggrqf(P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = INT( WORK( P+MN+1 ) );

      // Update c = Z**T *c = ( c1 ) N-P
                           // ( c2 ) M+P-N

      dormqr('Left', 'Transpose', M, 1, MN, A, LDA, WORK( P+1 ), C, max( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = max( LOPT, INT( WORK( P+MN+1 ) ) );

      // Solve T12*x2 = d for x2

      if ( P > 0 ) {
         dtrtrs('Upper', 'No transpose', 'Non-unit', P, 1, B( 1, N-P+1 ), LDB, D, P, INFO );

         if ( INFO > 0 ) {
            INFO = 1;
            return;
         }

         // Put the solution in X

         dcopy(P, D, 1, X( N-P+1 ), 1 );

         // Update c1

         dgemv('No transpose', N-P, P, -ONE, A( 1, N-P+1 ), LDA, D, 1, ONE, C, 1 );
      }

      // Solve R11*x1 = c1 for x1

      if ( N > P ) {
         dtrtrs('Upper', 'No transpose', 'Non-unit', N-P, 1, A, LDA, C, N-P, INFO );

         if ( INFO > 0 ) {
            INFO = 2;
            return;
         }

         // Put the solutions in X

         dcopy(N-P, C, 1, X, 1 );
      }

      // Compute the residual vector:

      if ( M < N ) {
         NR = M + P - N;
         if (NR > 0) dgemv( 'No transpose', NR, N-M, -ONE, A( N-P+1, M+1 ), LDA, D( NR+1 ), 1, ONE, C( N-P+1 ), 1 );
      } else {
         NR = P;
      }
      if ( NR > 0 ) {
         dtrmv('Upper', 'No transpose', 'Non unit', NR, A( N-P+1, N-P+1 ), LDA, D, 1 );
         daxpy(NR, -ONE, D, 1, C( N-P+1 ), 1 );
      }

      // Backward transformation x = Q**T*x

      dormrq('Left', 'Transpose', N, 1, P, B, LDB, WORK( 1 ), X, N, WORK( P+MN+1 ), LWORK-P-MN, INFO );
      WORK[1] = P + MN + max( LOPT, INT( WORK( P+MN+1 ) ) );

      return;
      }