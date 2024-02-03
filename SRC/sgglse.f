      SUBROUTINE SGGLSE( M, N, P, A, LDA, B, LDB, C, D, X, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), C( * ), D( * ), WORK( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKMIN, LWKOPT, MN, NB, NB1, NB2, NB3, NB4, NR;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGGRQF, SORMQR, SORMRQ, STRMV, STRTRS, XERBLA
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
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( P.LT.0 .OR. P.GT.N .OR. P.LT.N-M ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, P ) ) {
         INFO = -7
      }

      // Calculate workspace

      if ( INFO == 0) {
         if ( N == 0 ) {
            LWKMIN = 1
            LWKOPT = 1
         } else {
            NB1 = ILAENV( 1, 'SGEQRF', ' ', M, N, -1, -1 )
            NB2 = ILAENV( 1, 'SGERQF', ' ', M, N, -1, -1 )
            NB3 = ILAENV( 1, 'SORMQR', ' ', M, N, P, -1 )
            NB4 = ILAENV( 1, 'SORMRQ', ' ', M, N, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = P + MN + MAX( M, N )*NB
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

         if ( LWORK.LT.LWKMIN && .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO != 0 ) {
         xerbla('SGGLSE', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Compute the GRQ factorization of matrices B and A:

             // B*Q**T = (  0  T12 ) P   Z**T*A*Q**T = ( R11 R12 ) N-P
                         // N-P  P                     (  0  R22 ) M+P-N
                                                       // N-P  P

      // where T12 and R11 are upper triangular, and Q and Z are
      // orthogonal.

      sggrqf(P, M, N, B, LDB, WORK, A, LDA, WORK( P+1 ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = INT( WORK( P+MN+1 ) )

      // Update c = Z**T *c = ( c1 ) N-P
                           // ( c2 ) M+P-N

      sormqr('Left', 'Transpose', M, 1, MN, A, LDA, WORK( P+1 ), C, MAX( 1, M ), WORK( P+MN+1 ), LWORK-P-MN, INFO );
      LOPT = MAX( LOPT, INT( WORK( P+MN+1 ) ) )

      // Solve T12*x2 = d for x2

      if ( P.GT.0 ) {
         strtrs('Upper', 'No transpose', 'Non-unit', P, 1, B( 1, N-P+1 ), LDB, D, P, INFO );

         if ( INFO.GT.0 ) {
            INFO = 1
            RETURN
         }

         // Put the solution in X

         scopy(P, D, 1, X( N-P+1 ), 1 );

         // Update c1

         sgemv('No transpose', N-P, P, -ONE, A( 1, N-P+1 ), LDA, D, 1, ONE, C, 1 );
      }

      // Solve R11*x1 = c1 for x1

      if ( N.GT.P ) {
         strtrs('Upper', 'No transpose', 'Non-unit', N-P, 1, A, LDA, C, N-P, INFO );

         if ( INFO.GT.0 ) {
            INFO = 2
            RETURN
         }

         // Put the solutions in X

         scopy(N-P, C, 1, X, 1 );
      }

      // Compute the residual vector:

      if ( M.LT.N ) {
         NR = M + P - N
         if (NR.GT.0) CALL SGEMV( 'No transpose', NR, N-M, -ONE, A( N-P+1, M+1 ), LDA, D( NR+1 ), 1, ONE, C( N-P+1 ), 1 );
      } else {
         NR = P
      }
      if ( NR.GT.0 ) {
         strmv('Upper', 'No transpose', 'Non unit', NR, A( N-P+1, N-P+1 ), LDA, D, 1 );
         saxpy(NR, -ONE, D, 1, C( N-P+1 ), 1 );
      }

      // Backward transformation x = Q**T*x

      sormrq('Left', 'Transpose', N, 1, P, B, LDB, WORK( 1 ), X, N, WORK( P+MN+1 ), LWORK-P-MN, INFO );
      WORK( 1 ) = P + MN + MAX( LOPT, INT( WORK( P+MN+1 ) ) )

      RETURN

      // End of SGGLSE

      }
