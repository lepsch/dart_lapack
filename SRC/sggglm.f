      SUBROUTINE SGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), X( * ), Y( * )
      // ..

*  ===================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SGEMV, SGGQRF, SORMQR, SORMRQ, STRTRS, XERBLA
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
      NP = MIN( N, P )
      LQUERY = ( LWORK.EQ.-1 )
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( M.LT.0 .OR. M.GT.N ) {
         INFO = -2
      } else if ( P.LT.0 .OR. P.LT.N-M ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB.LT.MAX( 1, N ) ) {
         INFO = -7
      }

      // Calculate workspace

      if ( INFO.EQ.0) {
         if ( N.EQ.0 ) {
            LWKMIN = 1
            LWKOPT = 1
         } else {
            NB1 = ILAENV( 1, 'SGEQRF', ' ', N, M, -1, -1 )
            NB2 = ILAENV( 1, 'SGERQF', ' ', N, M, -1, -1 )
            NB3 = ILAENV( 1, 'SORMQR', ' ', N, M, P, -1 )
            NB4 = ILAENV( 1, 'SORMRQ', ' ', N, M, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = M + NP + MAX( N, P )*NB
         }
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT)

         if ( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SGGGLM', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N.EQ.0 ) {
         for (I = 1; I <= M; I++) {
            X(I) = ZERO
         END DO
         for (I = 1; I <= P; I++) {
            Y(I) = ZERO
         END DO
         RETURN
      }

      // Compute the GQR factorization of matrices A and B:

           // Q**T*A = ( R11 ) M,    Q**T*B*Z**T = ( T11   T12 ) M
                    // (  0  ) N-M                 (  0    T22 ) N-M
                       // M                         M+P-N  N-M

      // where R11 and T22 are upper triangular, and Q and Z are
      // orthogonal.

      sggqrf(N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = INT( WORK( M+NP+1 ) )

      // Update left-hand-side vector d = Q**T*d = ( d1 ) M
                                                // ( d2 ) N-M

      sormqr('Left', 'Transpose', N, 1, M, A, LDA, WORK, D, MAX( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = MAX( LOPT, INT( WORK( M+NP+1 ) ) )

      // Solve T22*y2 = d2 for y2

      if ( N.GT.M ) {
         strtrs('Upper', 'No transpose', 'Non unit', N-M, 1, B( M+1, M+P-N+1 ), LDB, D( M+1 ), N-M, INFO );

         if ( INFO.GT.0 ) {
            INFO = 1
            RETURN
         }

         scopy(N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 );
      }

      // Set y1 = 0

      DO 10 I = 1, M + P - N
         Y( I ) = ZERO
      } // 10

      // Update d1 = d1 - T12*y2

      sgemv('No transpose', M, N-M, -ONE, B( 1, M+P-N+1 ), LDB, Y( M+P-N+1 ), 1, ONE, D, 1 );

      // Solve triangular system: R11*x = d1

      if ( M.GT.0 ) {
         strtrs('Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D, M, INFO );

         if ( INFO.GT.0 ) {
            INFO = 2
            RETURN
         }

         // Copy D to X

         scopy(M, D, 1, X, 1 );
      }

      // Backward transformation y = Z**T *y

      sormrq('Left', 'Transpose', P, 1, NP, B( MAX( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y, MAX( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      WORK( 1 ) = M + NP + MAX( LOPT, INT( WORK( M+NP+1 ) ) )

      RETURN

      // End of SGGGLM

      }
