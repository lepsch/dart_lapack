      SUBROUTINE SSYGVD( ITYPE, JOBZ, UPLO, N, A, LDA, B, LDB, W, WORK, LWORK, IWORK, LIWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDA, LDB, LIWORK, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IWORK( * );
      REAL               A( LDA, * ), B( LDB, * ), W( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE
      const              ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER, WANTZ;
      String             TRANS;
      int                LIOPT, LIWMIN, LOPT, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL SPOTRF, SSYEVD, SSYGST, STRMM, STRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, REAL
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK == -1 || LIWORK == -1 )

      INFO = 0
      if ( N <= 1 ) {
         LIWMIN = 1
         LWMIN = 1
      } else if ( WANTZ ) {
         LIWMIN = 3 + 5*N
         LWMIN = 1 + 6*N + 2*N**2
      } else {
         LIWMIN = 1
         LWMIN = 2*N + 1
      }
      LOPT = LWMIN
      LIOPT = LIWMIN
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1
      } else if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( N < 0 ) {
         INFO = -4
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -6
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -8
      }

      if ( INFO == 0 ) {
         WORK( 1 ) = SROUNDUP_LWORK(LOPT)
         IWORK( 1 ) = LIOPT

         if ( LWORK < LWMIN && .NOT.LQUERY ) {
            INFO = -11
         } else if ( LIWORK < LIWMIN && .NOT.LQUERY ) {
            INFO = -13
         }
      }

      if ( INFO != 0 ) {
         xerbla('SSYGVD', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a Cholesky factorization of B.

      spotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      ssygst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      ssyevd(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, IWORK, LIWORK, INFO );
      LOPT = INT( MAX( REAL( LOPT ), REAL( WORK( 1 ) ) ) )
      LIOPT = INT( MAX( REAL( LIOPT ), REAL( IWORK( 1 ) ) ) )

      if ( WANTZ && INFO == 0 ) {

         // Backtransform eigenvectors to the original problem.

         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'T'
            }

            strsm('Left', UPLO, TRANS, 'Non-unit', N, N, ONE, B, LDB, A, LDA );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T'
            } else {
               TRANS = 'N'
            }

            strmm('Left', UPLO, TRANS, 'Non-unit', N, N, ONE, B, LDB, A, LDA );
         }
      }

      WORK( 1 ) = SROUNDUP_LWORK(LOPT)
      IWORK( 1 ) = LIOPT

      RETURN

      // End of SSYGVD

      }
