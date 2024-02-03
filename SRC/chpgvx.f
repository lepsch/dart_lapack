      SUBROUTINE CHPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO )

*  -- LAPACK driver routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             JOBZ, RANGE, UPLO;
      int                IL, INFO, ITYPE, IU, LDZ, M, N;
      REAL               ABSTOL, VL, VU
      // ..
      // .. Array Arguments ..
      int                IFAIL( * ), IWORK( * );
      REAL               RWORK( * ), W( * )
      COMPLEX            AP( * ), BP( * ), WORK( * ), Z( LDZ, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               ALLEIG, INDEIG, UPPER, VALEIG, WANTZ;
      String             TRANS;
      int                J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHPEVX, CHPGST, CPPTRF, CTPMV, CTPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      WANTZ = LSAME( JOBZ, 'V' )
      UPPER = LSAME( UPLO, 'U' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      INFO = 0
      if ( ITYPE.LT.1 || ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.( WANTZ || LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( ALLEIG || VALEIG || INDEIG ) ) {
         INFO = -3
      } else if ( .NOT.( UPPER || LSAME( UPLO, 'L' ) ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else {
         if ( VALEIG ) {
            if ( N.GT.0 && VU.LE.VL ) {
               INFO = -9
            }
         } else if ( INDEIG ) {
            if ( IL.LT.1 ) {
               INFO = -10
            } else if ( IU.LT.MIN( N, IL ) || IU.GT.N ) {
               INFO = -11
            }
         }
      }
      if ( INFO == 0 ) {
         if ( LDZ.LT.1 || ( WANTZ && LDZ.LT.N ) ) {
            INFO = -16
         }
      }

      if ( INFO != 0 ) {
         xerbla('CHPGVX', -INFO );
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
      chpevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, RWORK, IWORK, IFAIL, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         if (INFO.GT.0) M = INFO - 1;
         if ( ITYPE == 1 || ITYPE == 2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**H*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'C'
            }

            for (J = 1; J <= M; J++) { // 10
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

            for (J = 1; J <= M; J++) { // 20
               ctpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }

      RETURN

      // End of CHPGVX

      }
