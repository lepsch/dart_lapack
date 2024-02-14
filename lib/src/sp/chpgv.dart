      void chpgv(final int ITYPE, final int JOBZ, final int UPLO, final int N, final int AP, final int BP, final int W, final Matrix<double> Z_, final int LDZ, final Array<double> _WORK_, final Array<double> RWORK_, final Box<int> INFO,) {
  final Z = Z_.dim();
  final _WORK = _WORK_.dim();
  final RWORK = RWORK_.dim();

// -- LAPACK driver routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, N;
      double               RWORK( * ), W( * );
      Complex            AP( * ), BP( * ), WORK( * ), Z( LDZ, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      bool               UPPER, WANTZ;
      String             TRANS;
      int                J, NEIG;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPEV, CHPGST, CPPTRF, CTPMV, CTPSV, XERBLA

      // Test the input parameters.

      WANTZ = lsame( JOBZ, 'V' );
      UPPER = lsame( UPLO, 'U' );

      INFO = 0;
      if ( ITYPE < 1 || ITYPE > 3 ) {
         INFO = -1;
      } else if ( !( WANTZ || lsame( JOBZ, 'N' ) ) ) {
         INFO = -2;
      } else if ( !( UPPER || lsame( UPLO, 'L' ) ) ) {
         INFO = -3;
      } else if ( N < 0 ) {
         INFO = -4;
      } else if ( LDZ < 1 || ( WANTZ && LDZ < N ) ) {
         INFO = -9;
      }
      if ( INFO != 0 ) {
         xerbla('CHPGV ', -INFO );
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
      chpev(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N;
         if (INFO > 0) NEIG = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N';
            } else {
               TRANS = 'C';
            }

            for (J = 1; J <= NEIG; J++) { // 10
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

            for (J = 1; J <= NEIG; J++) { // 20
               ctpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }
      }
