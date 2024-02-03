      SUBROUTINE ZLAMTSQR( SIDE, TRANS, M, N, K, MB, NB, A, LDA, T, LDT, C, LDC, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE, TRANS;
      int                INFO, LDA, M, N, K, MB, NB, LDT, LWORK, LDC;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), WORK( * ), C( LDC, * ), T( LDT, * )
      // ..

* =====================================================================

      // ..
      // .. Local Scalars ..
      bool               LEFT, RIGHT, TRAN, NOTRAN, LQUERY;
      int                I, II, KK, LW, CTR, Q, MINMNK, LWMIN;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMQRT, ZTPMQRT, XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      LQUERY  = ( LWORK == -1 )
      NOTRAN  = LSAME( TRANS, 'N' )
      TRAN    = LSAME( TRANS, 'C' )
      LEFT    = LSAME( SIDE, 'L' )
      RIGHT   = LSAME( SIDE, 'R' )
      if ( LEFT ) {
        LW = N * NB
        Q = M
      } else {
        LW = M * NB
        Q = N
      }

      MINMNK = MIN( M, N, K )
      if ( MINMNK == 0 ) {
        LWMIN = 1
      } else {
        LWMIN = MAX( 1, LW )
      }

      if ( .NOT.LEFT .AND. .NOT.RIGHT ) {
        INFO = -1
      } else if ( .NOT.TRAN .AND. .NOT.NOTRAN ) {
        INFO = -2
      } else if ( M.LT.K ) {
        INFO = -3
      } else if ( N.LT.0 ) {
        INFO = -4
      } else if ( K.LT.0 ) {
        INFO = -5
      } else if ( K.LT.NB .OR. NB.LT.1 ) {
        INFO = -7
      } else if ( LDA.LT.MAX( 1, Q ) ) {
        INFO = -9
      } else if ( LDT.LT.MAX( 1, NB ) ) {
        INFO = -11
      } else if ( LDC.LT.MAX( 1, M ) ) {
        INFO = -13
      } else if ( LWORK.LT.LWMIN .AND. (.NOT.LQUERY) ) {
        INFO = -15
      }

      if ( INFO == 0 ) {
        WORK( 1 ) = LWMIN
      }

      if ( INFO.NE.0 ) {
        xerbla('ZLAMTSQR', -INFO );
        RETURN
      } else if ( LQUERY ) {
        RETURN
      }

      // Quick return if possible

      if ( MINMNK == 0 ) {
        RETURN
      }

      // Determine the block size if it is tall skinny or short and wide

      if ((MB.LE.K).OR.(MB.GE.MAX(M,N,K))) {
        zgemqrt(SIDE, TRANS, M, N, K, NB, A, LDA, T, LDT, C, LDC, WORK, INFO );
        RETURN
      }

      if (LEFT.AND.NOTRAN) {

          // Multiply Q to the last block of C

         KK = MOD((M-K),(MB-K))
         CTR = (M-K)/(MB-K)
         if (KK.GT.0) {
           II=M-KK+1
           ztpmqrt('L','N',KK , N, K, 0, NB, A(II,1), LDA, T(1, CTR * K + 1),LDT , C(1,1), LDC, C(II,1), LDC, WORK, INFO );
         } else {
           II=M+1
         }

         DO I=II-(MB-K),MB+1,-(MB-K)

          // Multiply Q to the current block of C (I:I+MB,1:N)

           CTR = CTR - 1
           ztpmqrt('L','N',MB-K , N, K, 0,NB, A(I,1), LDA, T(1,CTR * K + 1),LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO );

         }

          // Multiply Q to the first block of C (1:MB,1:N)

         zgemqrt('L','N',MB , N, K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

      } else if (LEFT.AND.TRAN) {

          // Multiply Q to the first block of C

         KK = MOD((M-K),(MB-K))
         II=M-KK+1
         CTR = 1
         zgemqrt('L','C',MB , N, K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

         DO I=MB+1,II-MB+K,(MB-K)

          // Multiply Q to the current block of C (I:I+MB,1:N)

          ztpmqrt('L','C',MB-K , N, K, 0,NB, A(I,1), LDA, T(1,CTR * K + 1),LDT, C(1,1), LDC, C(I,1), LDC, WORK, INFO );
          CTR = CTR + 1

         }
         if (II.LE.M) {

          // Multiply Q to the last block of C

          ztpmqrt('L','C',KK , N, K, 0,NB, A(II,1), LDA, T(1, CTR * K + 1), LDT, C(1,1), LDC, C(II,1), LDC, WORK, INFO );

         }

      } else if (RIGHT.AND.TRAN) {

          // Multiply Q to the last block of C

          KK = MOD((N-K),(MB-K))
          CTR = (N-K)/(MB-K)
          if (KK.GT.0) {
            II=N-KK+1
            ztpmqrt('R','C',M , KK, K, 0, NB, A(II,1), LDA, T(1,CTR * K + 1), LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO );
          } else {
            II=N+1
          }

          DO I=II-(MB-K),MB+1,-(MB-K)

          // Multiply Q to the current block of C (1:M,I:I+MB)

            CTR = CTR - 1
            ztpmqrt('R','C',M , MB-K, K, 0,NB, A(I,1), LDA, T(1, CTR * K + 1), LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO );

          }

          // Multiply Q to the first block of C (1:M,1:MB)

          zgemqrt('R','C',M , MB, K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

      } else if (RIGHT.AND.NOTRAN) {

          // Multiply Q to the first block of C

         KK = MOD((N-K),(MB-K))
         II=N-KK+1
         CTR = 1
         zgemqrt('R','N', M, MB , K, NB, A(1,1), LDA, T ,LDT ,C(1,1), LDC, WORK, INFO );

         DO I=MB+1,II-MB+K,(MB-K)

          // Multiply Q to the current block of C (1:M,I:I+MB)

          ztpmqrt('R','N', M, MB-K, K, 0,NB, A(I,1), LDA, T(1, CTR * K + 1),LDT, C(1,1), LDC, C(1,I), LDC, WORK, INFO );
          CTR = CTR + 1

         }
         if (II.LE.N) {

          // Multiply Q to the last block of C

          ztpmqrt('R','N', M, KK , K, 0,NB, A(II,1), LDA, T(1,CTR * K + 1),LDT, C(1,1), LDC, C(1,II), LDC, WORK, INFO );

         }

      }

      WORK( 1 ) = LWMIN
      RETURN

      // End of ZLAMTSQR

      }
