      SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TYPE;
      int                INFO, KL, KU, LDA, M, N;
      double             CFROM, CTO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = 0.0D0, ONE = 1.0D0 ;
      // ..
      // .. Local Scalars ..
      bool               DONE;
      int                I, ITYPE, J, K1, K2, K3, K4;
      double             BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      double             DLAMCH;
      // EXTERNAL LSAME, DLAMCH, DISNAN
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0

      if ( LSAME( TYPE, 'G' ) ) {
         ITYPE = 0
      } else if ( LSAME( TYPE, 'L' ) ) {
         ITYPE = 1
      } else if ( LSAME( TYPE, 'U' ) ) {
         ITYPE = 2
      } else if ( LSAME( TYPE, 'H' ) ) {
         ITYPE = 3
      } else if ( LSAME( TYPE, 'B' ) ) {
         ITYPE = 4
      } else if ( LSAME( TYPE, 'Q' ) ) {
         ITYPE = 5
      } else if ( LSAME( TYPE, 'Z' ) ) {
         ITYPE = 6
      } else {
         ITYPE = -1
      }

      if ( ITYPE.EQ.-1 ) {
         INFO = -1
      } else if ( CFROM.EQ.ZERO .OR. DISNAN(CFROM) ) {
         INFO = -4
      } else if ( DISNAN(CTO) ) {
         INFO = -5
      } else if ( M.LT.0 ) {
         INFO = -6
      } else if ( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR. ( ITYPE.EQ.5 .AND. N.NE.M ) ) {
         INFO = -7
      } else if ( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) {
         INFO = -9
      } else if ( ITYPE.GE.4 ) {
         if ( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) {
            INFO = -2
         } else if ( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR. ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) ) {
            INFO = -3
         } else if ( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR. ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR. ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) {
            INFO = -9
         }
      }

      if ( INFO.NE.0 ) {
         xerbla('ZLASCL', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N.EQ.0 .OR. M.EQ.0) RETURN;

      // Get machine parameters

      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

      CFROMC = CFROM
      CTOC = CTO

      } // 10
      CFROM1 = CFROMC*SMLNUM
      if ( CFROM1.EQ.CFROMC ) {
         // CFROMC is an inf.  Multiply by a correctly signed zero for
         // finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = true;
         CTO1 = CTOC
      } else {
         CTO1 = CTOC / BIGNUM
         if ( CTO1.EQ.CTOC ) {
            // CTOC is either 0 or an inf.  In both cases, CTOC itself
            // serves as the correct multiplication factor.
            MUL = CTOC
            DONE = true;
            CFROMC = ONE
         } else if ( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) {
            MUL = SMLNUM
            DONE = false;
            CFROMC = CFROM1
         } else if ( ABS( CTO1 ).GT.ABS( CFROMC ) ) {
            MUL = BIGNUM
            DONE = false;
            CTOC = CTO1
         } else {
            MUL = CTOC / CFROMC
            DONE = true;
            if (MUL .EQ. ONE) RETURN;
         }
      }

      if ( ITYPE.EQ.0 ) {

         // Full matrix

         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= M; I++) { // 20
               A( I, J ) = A( I, J )*MUL
            } // 20
         } // 30

      } else if ( ITYPE.EQ.1 ) {

         // Lower triangular matrix

         for (J = 1; J <= N; J++) { // 50
            for (I = J; I <= M; I++) { // 40
               A( I, J ) = A( I, J )*MUL
            } // 40
         } // 50

      } else if ( ITYPE.EQ.2 ) {

         // Upper triangular matrix

         for (J = 1; J <= N; J++) { // 70
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
            } // 60
         } // 70

      } else if ( ITYPE.EQ.3 ) {

         // Upper Hessenberg matrix

         for (J = 1; J <= N; J++) { // 90
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
            } // 80
         } // 90

      } else if ( ITYPE.EQ.4 ) {

         // Lower half of a symmetric band matrix

         K3 = KL + 1
         K4 = N + 1
         for (J = 1; J <= N; J++) { // 110
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
            } // 100
         } // 110

      } else if ( ITYPE.EQ.5 ) {

         // Upper half of a symmetric band matrix

         K1 = KU + 2
         K3 = KU + 1
         for (J = 1; J <= N; J++) { // 130
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
            } // 120
         } // 130

      } else if ( ITYPE.EQ.6 ) {

         // Band matrix

         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         for (J = 1; J <= N; J++) { // 150
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
            } // 140
         } // 150

      }

      if (.NOT.DONE) GO TO 10;

      RETURN

      // End of ZLASCL

      }
