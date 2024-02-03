      SUBROUTINE CLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TYPE;
      int                INFO, KL, KU, LDA, M, N;
      REAL               CFROM, CTO
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E0, ONE = 1.0E0 ;
      // ..
      // .. Local Scalars ..
      bool               DONE;
      int                I, ITYPE, J, K1, K2, K3, K4;
      REAL               BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
      // ..
      // .. External Functions ..
      bool               LSAME, SISNAN;
      REAL               SLAMCH
      // EXTERNAL LSAME, SLAMCH, SISNAN
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

      if ( ITYPE == -1 ) {
         INFO = -1
      } else if ( CFROM == ZERO || SISNAN(CFROM) ) {
         INFO = -4
      } else if ( SISNAN(CTO) ) {
         INFO = -5
      } else if ( M < 0 ) {
         INFO = -6
      } else if ( N < 0 || ( ITYPE == 4 && N != M ) || ( ITYPE == 5 && N != M ) ) {
         INFO = -7
      } else if ( ITYPE.LE.3 && LDA < MAX( 1, M ) ) {
         INFO = -9
      } else if ( ITYPE >= 4 ) {
         if ( KL < 0 || KL > MAX( M-1, 0 ) ) {
            INFO = -2
         } else if ( KU < 0 || KU > MAX( N-1, 0 ) || ( ( ITYPE == 4 || ITYPE == 5 ) && KL != KU ) ) {
            INFO = -3
         } else if ( ( ITYPE == 4 && LDA < KL+1 ) || ( ITYPE == 5 && LDA < KU+1 ) || ( ITYPE == 6 && LDA < 2*KL+KU+1 ) ) {
            INFO = -9
         }
      }

      if ( INFO != 0 ) {
         xerbla('CLASCL', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0 || M == 0) RETURN;

      // Get machine parameters

      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

      CFROMC = CFROM
      CTOC = CTO

      } // 10
      CFROM1 = CFROMC*SMLNUM
      if ( CFROM1 == CFROMC ) {
         // CFROMC is an inf.  Multiply by a correctly signed zero for
         // finite CTOC, or a NaN if CTOC is infinite.
         MUL = CTOC / CFROMC
         DONE = true;
         CTO1 = CTOC
      } else {
         CTO1 = CTOC / BIGNUM
         if ( CTO1 == CTOC ) {
            // CTOC is either 0 or an inf.  In both cases, CTOC itself
            // serves as the correct multiplication factor.
            MUL = CTOC
            DONE = true;
            CFROMC = ONE
         } else if ( ABS( CFROM1 ) > ABS( CTOC ) && CTOC != ZERO ) {
            MUL = SMLNUM
            DONE = false;
            CFROMC = CFROM1
         } else if ( ABS( CTO1 ) > ABS( CFROMC ) ) {
            MUL = BIGNUM
            DONE = false;
            CTOC = CTO1
         } else {
            MUL = CTOC / CFROMC
            DONE = true;
            if (MUL == ONE) RETURN;
         }
      }

      if ( ITYPE == 0 ) {

         // Full matrix

         for (J = 1; J <= N; J++) { // 30
            for (I = 1; I <= M; I++) { // 20
               A( I, J ) = A( I, J )*MUL
            } // 20
         } // 30

      } else if ( ITYPE == 1 ) {

         // Lower triangular matrix

         for (J = 1; J <= N; J++) { // 50
            for (I = J; I <= M; I++) { // 40
               A( I, J ) = A( I, J )*MUL
            } // 40
         } // 50

      } else if ( ITYPE == 2 ) {

         // Upper triangular matrix

         for (J = 1; J <= N; J++) { // 70
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
            } // 60
         } // 70

      } else if ( ITYPE == 3 ) {

         // Upper Hessenberg matrix

         for (J = 1; J <= N; J++) { // 90
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
            } // 80
         } // 90

      } else if ( ITYPE == 4 ) {

         // Lower half of a symmetric band matrix

         K3 = KL + 1
         K4 = N + 1
         for (J = 1; J <= N; J++) { // 110
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
            } // 100
         } // 110

      } else if ( ITYPE == 5 ) {

         // Upper half of a symmetric band matrix

         K1 = KU + 2
         K3 = KU + 1
         for (J = 1; J <= N; J++) { // 130
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
            } // 120
         } // 130

      } else if ( ITYPE == 6 ) {

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

      // End of CLASCL

      }
