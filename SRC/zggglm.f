      SUBROUTINE ZGGGLM( N, M, P, A, LDA, B, LDB, D, X, Y, WORK, LWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), D( * ), WORK( * ), X( * ), Y( * )
      // ..

*  ===================================================================

      // .. Parameters ..
      COMPLEX*16         CZERO, CONE
      const              CZERO = ( 0.0D+0, 0.0D+0 ), CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, LOPT, LWKMIN, LWKOPT, NB, NB1, NB2, NB3, NB4, NP;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZCOPY, ZGEMV, ZGGQRF, ZTRTRS, ZUNMQR, ZUNMRQ
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      NP = MIN( N, P )
      LQUERY = ( LWORK == -1 )
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

      if ( INFO == 0) {
         if ( N == 0 ) {
            LWKMIN = 1
            LWKOPT = 1
         } else {
            NB1 = ILAENV( 1, 'ZGEQRF', ' ', N, M, -1, -1 )
            NB2 = ILAENV( 1, 'ZGERQF', ' ', N, M, -1, -1 )
            NB3 = ILAENV( 1, 'ZUNMQR', ' ', N, M, P, -1 )
            NB4 = ILAENV( 1, 'ZUNMRQ', ' ', N, M, P, -1 )
            NB = MAX( NB1, NB2, NB3, NB4 )
            LWKMIN = M + N + P
            LWKOPT = M + NP + MAX( N, P )*NB
         }
         WORK( 1 ) = LWKOPT

         if ( LWORK.LT.LWKMIN && .NOT.LQUERY ) {
            INFO = -12
         }
      }

      if ( INFO != 0 ) {
         xerbla('ZGGGLM', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if ( N == 0 ) {
         for (I = 1; I <= M; I++) {
            X(I) = CZERO
         }
         for (I = 1; I <= P; I++) {
            Y(I) = CZERO
         }
         RETURN
      }

      // Compute the GQR factorization of matrices A and B:

           // Q**H*A = ( R11 ) M,    Q**H*B*Z**H = ( T11   T12 ) M
                    // (  0  ) N-M                 (  0    T22 ) N-M
                       // M                         M+P-N  N-M

      // where R11 and T22 are upper triangular, and Q and Z are
      // unitary.

      zggqrf(N, M, P, A, LDA, WORK, B, LDB, WORK( M+1 ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = INT( WORK( M+NP+1 ) )

      // Update left-hand-side vector d = Q**H*d = ( d1 ) M
                                                // ( d2 ) N-M

      zunmqr('Left', 'Conjugate transpose', N, 1, M, A, LDA, WORK, D, MAX( 1, N ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      LOPT = MAX( LOPT, INT( WORK( M+NP+1 ) ) )

      // Solve T22*y2 = d2 for y2

      if ( N.GT.M ) {
         ztrtrs('Upper', 'No transpose', 'Non unit', N-M, 1, B( M+1, M+P-N+1 ), LDB, D( M+1 ), N-M, INFO );

         if ( INFO.GT.0 ) {
            INFO = 1
            RETURN
         }

         zcopy(N-M, D( M+1 ), 1, Y( M+P-N+1 ), 1 );
      }

      // Set y1 = 0

      for (I = 1; I <= M + P - N; I++) { // 10
         Y( I ) = CZERO
      } // 10

      // Update d1 = d1 - T12*y2

      zgemv('No transpose', M, N-M, -CONE, B( 1, M+P-N+1 ), LDB, Y( M+P-N+1 ), 1, CONE, D, 1 );

      // Solve triangular system: R11*x = d1

      if ( M.GT.0 ) {
         ztrtrs('Upper', 'No Transpose', 'Non unit', M, 1, A, LDA, D, M, INFO );

         if ( INFO.GT.0 ) {
            INFO = 2
            RETURN
         }

         // Copy D to X

         zcopy(M, D, 1, X, 1 );
      }

      // Backward transformation y = Z**H *y

      zunmrq('Left', 'Conjugate transpose', P, 1, NP, B( MAX( 1, N-P+1 ), 1 ), LDB, WORK( M+1 ), Y, MAX( 1, P ), WORK( M+NP+1 ), LWORK-M-NP, INFO );
      WORK( 1 ) = M + NP + MAX( LOPT, INT( WORK( M+NP+1 ) ) )

      RETURN

      // End of ZGGGLM

      }
