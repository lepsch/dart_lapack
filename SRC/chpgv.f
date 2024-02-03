      SUBROUTINE CHPGV( ITYPE, JOBZ, UPLO, N, AP, BP, W, Z, LDZ, WORK, RWORK, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, UPLO;
      int                INFO, ITYPE, LDZ, N;
      // ..
      // .. Array Arguments ..
      REAL               RWORK( * ), W( * )
      COMPLEX            AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               UPPER, WANTZ;
      String             TRANS;
      int                J, NEIG;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPEV, CHPGST, CPPTRF, CTPMV, CTPSV, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )

      INFO = 0
      if ( ITYPE.LT.1 || ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -3
      } else if ( N.LT.0 ) {
         INFO = -4
      } else if ( LDZ.LT.1 || ( WANTZ && LDZ.LT.N ) ) {
         INFO = -9
      }
      if ( INFO != 0 ) {
         xerbla('CHPGV ', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Form a Cholesky factorization of B.

      cpptrf(UPLO, N, BP, INFO );
      if ( INFO != 0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      chpgst(ITYPE, UPLO, N, AP, BP, INFO );
      chpev(JOBZ, UPLO, N, AP, W, Z, LDZ, WORK, RWORK, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         NEIG = N
         if (INFO.GT.0) NEIG = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'C'
            }

            for (J = 1; J <= NEIG; J++) { // 10
               ctpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 10

         } else if ( ITYPE == 3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**H*y

            if ( UPPER ) {
               TRANS = 'C'
            } else {
               TRANS = 'N'
            }

            for (J = 1; J <= NEIG; J++) { // 20
               ctpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }
      RETURN

      // End of CHPGV

      }
