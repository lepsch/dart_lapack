      SUBROUTINE CGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( * ), D( * ), WORK( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4, NR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CGGRQF, CTRMV, CTRTRS, CUNMQR, CUNMRQ, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      MN = MIN( M, N )
      LQUERY = ( LWORK == -1 )
      if ( M < 0 ) {
         INFO = -1
      } else if ( N < 0 ) {
         INFO = -2
      } else if ( P < 0 || P > N || P < N-M ) {
         INFO = -3
      } else if ( LDA < MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, P ) ) {
         INFO = -7
      }

      // Calculate workspace

      if ( INFO == 0) {
         if ( N == 0 ) {
            LWKMIN = 1
            LWKOPT = 1
         } else {
            NB1 = ILAENV( 1, 'CGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'CGERQF', ' ', M, N, -1, -1 )
            NB3 = ILAENV( 1, 'CUNMQR', ' ', M, N, P, -1 )
            NB4 = ILAENV( 1, 'CUNMRQ', ' ', M, N, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = P + MN + MAX( M, N )*NB
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

         if ( LWORK < LWKMIN && .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO != 0 ) {
         xerbla('CGGLSE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Compute the GRQ factorization of matrices B and A:

             // B*Q**H = (  0  T12 ) P   Z**H*A*Q**H = ( R11 R12 ) N-P
                         // N-P  P                     (  0  R22 ) M+P-N
                                                       // N-P  P

      // where T12 and R11 are upper triangular, and Q and Z are
      // unitary.

      cggrqf(P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = INT( WORK( P+MN+1 ) )

      // Update c = Z**H *c = ( c1 ) N-P
                        // ( c2 ) M+P-N

      cunmqr('Left', 'Conjugate Transpose', M, 1, MN, A, LDA, WORK( P+1 ), C, MAX( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = MAX( LOPT, INT( WORK( P+MN+1 ) ) )

      // Solve T12*x2 = d for x2

      if ( P > 0 ) {
         ctrtrs('Upper', 'No transpose', 'Non-unit', P, 1, B( 1, N-P+1 ), LDB, D, P, INFO );

         if ( INFO > 0 ) {
            INFO = 1
            RETURN
         }

         // Put the solution in X

      ccopy(P, D, 1, X( N-P+1 ), 1 );

         // Update c1

         cgemv('No transpose', N-P, P, -CONE, A( 1, N-P+1 ), LDA, D, 1, CONE, C, 1 );
      }

      // Solve R11*x1 = c1 for x1

      if ( N > P ) {
         ctrtrs('Upper', 'No transpose', 'Non-unit', N-P, 1, A, LDA, C, N-P, INFO );

         if ( INFO > 0 ) {
            INFO = 2
            RETURN
         }

         // Put the solutions in X

         ccopy(N-P, C, 1, X, 1 );
      }

      // Compute the residual vector:

      if ( M < N ) {
         NR = M + P - N
         if (NR > 0) CALL CGEMV( 'No transpose', NR, N-M, -CONE, A( N-P+1, M+1 ), LDA, D( NR+1 ), 1, CONE, C( N-P+1 ), 1 );
      } else {
         NR = P
      }
      if ( NR > 0 ) {
         ctrmv('Upper', 'No transpose', 'Non unit', NR, A( N-P+1, N-P+1 ), LDA, D, 1 );
         caxpy(NR, -CONE, D, 1, C( N-P+1 ), 1 );
      }

      // Backward transformation x = Q**H*x

      cunmrq('Left', 'Conjugate Transpose', N, 1, P, B, LDB, WORK( 1 ), X, N, WORK( P+MN+1 ), LWORK-P-MN, INFO );
      WORK( 1 ) = P + MN + MAX( LOPT, INT( WORK( P+MN+1 ) ) )

      RETURN

      // End of CGGLSE

      }
