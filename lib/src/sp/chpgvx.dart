      void chpgvx(ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO ) {

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, ITYPE, IU, LDZ, M, N;
      REAL               ABSTOL, VL, VU;
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               RWORK( * ), W( * );
      Complex            AP( * ), BP( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, UPPER, VALEIG, WANTZ;
      String             TRANS;
      int                J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPEVX, CHPGST, CPPTRF, CTPMV, CTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );
      ALLEIG = lsame( RANGE, 'A' );
      VALEIG = lsame( RANGE, 'V' );
      INDEIG = lsame( RANGE, 'I' );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -3;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -4;
      } else if ( N < 0 ) {
         INFO = -5;
      } else {
         if ( VALEIG ) {
            if ( N > 0 && VU <= VL ) {
               INFO = -9;
            }
         } else if ( INDEIG ) {
            if ( IL < 1 ) {
               INFO = -10;
            } else if ( IU < min( N, IL ) || IU > N ) {
               INFO = -11;
            }
         }
      }
      if ( INFO == 0 ) {
         if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
            INFO = -16;
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHPGVX', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Form a Cholesky factorization of B.

      cpptrf(UPLO, N, BP, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO;
         return;
      }

      // Transform problem to standard eigenvalue problem and solve.

      chpgst(ITYPE, UPLO, N, AP, BP, INFO );
      chpevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO );

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

            for (J = 1; J <= M; J++) { // 10
               ctpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 10

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H*y

            if ( UPPER ) {
               TRANS = 'C';
            } else {
               TRANS = 'N';
            }

            for (J = 1; J <= M; J++) { // 20
               ctpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }

      return;
      }
