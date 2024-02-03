      SUBROUTINE CHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO );

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N;
      REAL               ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               RWORK( * ), W( * );
      COMPLEX            A( LDA, * ), B( LDB, * ), WORK( * ), Z( LDZ, * );
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, LQUERY, UPPER, VALEIG, WANTZ;
      String             TRANS;
      int                LWKOPT, NB;
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV;
      REAL               SROUNDUP_LWORK;
      // EXTERNAL ILAENV, LSAME, SROUNDUP_LWORK
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHEEVX, CHEGST, CPOTRF, CTRMM, CTRSM, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' );
      UPPER = LSAME( UPLO, 'U' );
      ALLEIG = LSAME( RANGE, 'A' );
      VALEIG = LSAME( RANGE, 'V' );
      INDEIG = LSAME( RANGE, 'I' );
      LQUERY = ( LWORK == -1 );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -3;
      } else if ( !( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else if ( LDA < MAX( 1, N ) ) {
         INFO = -7;
      } else if ( LDB < MAX( 1, N ) ) {
         INFO = -9;
      } else {
         if ( VALEIG ) {
            if (N > 0 && VU <= VL) INFO = -11;
         } else if ( INDEIG ) {
            if ( IL < 1 || IL > MAX( 1, N ) ) {
               INFO = -12;
            } else if ( IU < MIN( N, IL ) || IU > N ) {
               INFO = -13;
            }
         }
      }
      if (INFO == 0) {
         if (LDZ < 1 || (WANTZ && LDZ < N)) {
            INFO = -18;
         }
      }

      if ( INFO == 0 ) {
         NB = ILAENV( 1, 'CHETRD', UPLO, N, -1, -1, -1 );
         LWKOPT = MAX( 1, ( NB + 1 )*N );
         WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);

         if ( LWORK < MAX( 1, 2*N ) && !LQUERY ) {
            INFO = -20;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHEGVX', -INFO );
         RETURN;
      } else if ( LQUERY ) {
         RETURN;
      }

      // Quick return if possible

      M = 0;
      if ( N == 0 ) {
         RETURN;
      }

      // Form a Cholesky factorization of B.

      cpotrf(UPLO, N, B, LDB, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         RETURN;
      }

      // Transform problem to standard eigenvalue problem and solve.

      chegst(ITYPE, UPLO, N, A, LDA, B, LDB, INFO );
      cheevx(JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, IFAIL, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         if (INFO > 0) M = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            ctrsm('Left', UPLO, TRANS, 'Non-unit', N, M, CONE, B, LDB, Z, LDZ );

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H*y

            if ( UPPER ) {
               TRANS = 'C';
            } else {
               TRANS = 'N';
            }

            ctrmm('Left', UPLO, TRANS, 'Non-unit', N, M, CONE, B, LDB, Z, LDZ );
         }
      }

      // Set WORK(1) to optimal complex workspace size.

      WORK( 1 ) = SROUNDUP_LWORK(LWKOPT);

      RETURN;

      // End of CHEGVX

      }
