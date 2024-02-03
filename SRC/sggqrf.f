      SUBROUTINE SGGQRF( N, M, P, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      REAL               A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQRF, SGERQF, SORMQR, XERBLA
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      REAL               SROUNDUP_LWORK
      // EXTERNAL SROUNDUP_LWORK
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      NB1 = ILAENV( 1, 'SGEQRF', ' ', N, M, -1, -1 )
      NB2 = ILAENV( 1, 'SGERQF', ' ', N, P, -1, -1 )
      NB3 = ILAENV( 1, 'SORMQR', ' ', N, M, P, -1 )
      NB = MAX( NB1, NB2, NB3 )
      LWKOPT = MAX( 1, MAX( N, M, P )*NB )
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      LQUERY = ( LWORK == -1 )
      if ( N < 0 ) {
         INFO = -1
      } else if ( M < 0 ) {
         INFO = -2
      } else if ( P < 0 ) {
         INFO = -3
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8
      } else if ( LWORK < MAX( 1, N, M, P ) && !LQUERY ) {
         INFO = -11
      }
      if ( INFO != 0 ) {
         xerbla('SGGQRF', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // QR factorization of N-by-M matrix A: A = Q*R

      sgeqrf(N, M, A, LDA, TAUA, WORK, LWORK, INFO );
      LOPT = INT( WORK( 1 ) )

      // Update B := Q**T*B.

      sormqr('Left', 'Transpose', N, P, MIN( N, M ), A, LDA, TAUA, B, LDB, WORK, LWORK, INFO );
      LOPT = MAX( LOPT, INT( WORK( 1 ) ) )

      // RQ factorization of N-by-P matrix B: B = T*Z.

      sgerqf(N, P, B, LDB, TAUB, WORK, LWORK, INFO );
      LWKOPT = MAX( LOPT, INT( WORK( 1 ) ) )

      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )

      RETURN

      // End of SGGQRF

      }
