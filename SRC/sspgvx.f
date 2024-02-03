      SUBROUTINE SSPGVX( ITYPE, JOBZ, RANGE, UPLO, N, AP, BP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO )

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
      REAL               AP( * ), BP( * ), W( * ), WORK( * ), Z( LDZ, * )
      // ..

* =====================================================================

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
      // EXTERNAL SPPTRF, SSPEVX, SSPGST, STPMV, STPSV, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      UPPER = LSAME( UPLO, 'U' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )

      INFO = 0
      if ( ITYPE.LT.1 .OR. ITYPE.GT.3 ) {
         INFO = -1
      } else if ( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) {
         INFO = -2
      } else if ( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) {
         INFO = -3
      } else if ( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -4
      } else if ( N.LT.0 ) {
         INFO = -5
      } else {
         if ( VALEIG ) {
            if ( N.GT.0 .AND. VU.LE.VL ) {
               INFO = -9
            }
         } else if ( INDEIG ) {
            if ( IL.LT.1 ) {
               INFO = -10
            } else if ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) {
               INFO = -11
            }
         }
      }
      if ( INFO.EQ.0 ) {
         if ( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) {
            INFO = -16
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('SSPGVX', -INFO );
         RETURN
      }

      // Quick return if possible

      M = 0
      IF( N.EQ.0 ) RETURN

      // Form a Cholesky factorization of B.

      spptrf(UPLO, N, BP, INFO );
      if ( INFO.NE.0 ) {
         INFO = N + INFO
         RETURN
      }

      // Transform problem to standard eigenvalue problem and solve.

      sspgst(ITYPE, UPLO, N, AP, BP, INFO );
      sspevx(JOBZ, RANGE, UPLO, N, AP, VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK, IWORK, IFAIL, INFO );

      if ( WANTZ ) {

         // Backtransform eigenvectors to the original problem.

         IF( INFO.GT.0 ) M = INFO - 1
         if ( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) {

            // For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
            // backtransform eigenvectors: x = inv(L)**T*y or inv(U)*y

            if ( UPPER ) {
               TRANS = 'N'
            } else {
               TRANS = 'T'
            }

            for (J = 1; J <= M; J++) { // 10
               stpsv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 10

         } else if ( ITYPE.EQ.3 ) {

            // For B*A*x=(lambda)*x;
            // backtransform eigenvectors: x = L*y or U**T*y

            if ( UPPER ) {
               TRANS = 'T'
            } else {
               TRANS = 'N'
            }

            for (J = 1; J <= M; J++) { // 20
               stpmv(UPLO, TRANS, 'Non-unit', N, BP, Z( 1, J ), 1 );
            } // 20
         }
      }

      RETURN

      // End of SSPGVX

      }
