      void zhetrd_he2hb(final int UPLO, final int N, final int KD, final Matrix<double> A, final int LDA, final Matrix<double> AB, final int LDAB, final int TAU, final Array<double> WORK, final int LWORK, final Box<int> INFO,) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      String             UPLO;
      int                INFO, LDA, LDAB, LWORK, N, KD;
      Complex         A( LDA, * ), AB( LDAB, * ),  TAU( * ), WORK( * );
      // ..

      double             RONE;
      Complex         ZERO, ONE, HALF;
      const              RONE = 1.0, ZERO = ( 0.0, 0.0 ), ONE = ( 1.0, 0.0 ), HALF = ( 0.5, 0.0 ) ;
      bool               LQUERY, UPPER;
      int                I, J, IINFO, LWMIN, PN, PK, LK, LDT, LDW, LDS2, LDS1, LS2, LS1, LW, LT, TPOS, WPOS, S2POS, S1POS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZHER2K, ZHEMM, ZGEMM, ZCOPY, ZLARFT, ZGELQF, ZGEQRF, ZLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- int                ILAENV2STAGE;
      // EXTERNAL lsame, ILAENV2STAGE

      // Determine the minimal workspace size required
      // and test the input parameters

      INFO   = 0;
      UPPER  = lsame( UPLO, 'U' );
      LQUERY = ( LWORK == -1 );
      if ( N <= KD+1 ) {
         LWMIN = 1;
      } else {
         LWMIN = ILAENV2STAGE( 4, 'ZHETRD_HE2HB', '', N, KD, -1, -1 );
      }

      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( KD < 0 ) {
         INFO = -3;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -5;
      } else if ( LDAB < max( 1, KD+1 ) ) {
         INFO = -7;
      } else if ( LWORK < LWMIN && !LQUERY ) {
         INFO = -10;
      }

      if ( INFO != 0 ) {
         xerbla('ZHETRD_HE2HB', -INFO );
         return;
      } else if ( LQUERY ) {
         WORK[1] = LWMIN;
         return;
      }

      // Quick return if possible
      // Copy the upper/lower portion of A into AB

      if ( N <= KD+1 ) {
          if ( UPPER ) {
              for (I = 1; I <= N; I++) { // 100
                  LK = min( KD+1, I );
                  zcopy(LK, A( I-LK+1, I ), 1,  AB( KD+1-LK+1, I ), 1 );
              } // 100
          } else {
              for (I = 1; I <= N; I++) { // 110
                  LK = min( KD+1, N-I+1 );
                  zcopy(LK, A( I, I ), 1, AB( 1, I ), 1 );
              } // 110
          }
          WORK[1] = 1;
          return;
      }

      // Determine the pointer position for the workspace

      LDT    = KD;
      LDS1   = KD;
      LT     = LDT*KD;
      LW     = N*KD;
      LS1    = LDS1*KD;
      LS2    = LWMIN - LT - LW - LS1;
       // LS2 = N*max(KD,FACTOPTNB)
      TPOS   = 1;
      WPOS   = TPOS  + LT;
      S1POS  = WPOS  + LW;
      S2POS  = S1POS + LS1;
      if ( UPPER ) {
          LDW    = KD;
          LDS2   = KD;
      } else {
          LDW    = N;
          LDS2   = N;
      }


      // Set the workspace of the triangular matrix T to zero once such a
      // way every time T is generated the upper/lower portion will be always zero

      zlaset("A", LDT, KD, ZERO, ZERO, WORK( TPOS ), LDT );

      if ( UPPER ) {
          for (I = 1; KD < 0 ? I >= N - KD : I <= N - KD; I += KD) { // 10
             PN = N-I-KD+1;
             PK = min( N-I-KD+1, KD );

             // Compute the LQ factorization of the current block

             zgelqf(KD, PN, A( I, I+KD ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO );

             // Copy the upper portion of A into AB

             for (J = I; J <= I+PK-1; J++) { // 20
                LK = min( KD, N-J ) + 1;
                zcopy(LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 );
             } // 20

             zlaset('Lower', PK, PK, ZERO, ONE,  A( I, I+KD ), LDA );

             // Form the matrix T

             zlarft('Forward', 'Rowwise', PN, PK, A( I, I+KD ), LDA, TAU( I ), WORK( TPOS ), LDT );

             // Compute W:

             zgemm('Conjugate', 'No transpose', PK, PN, PK, ONE,  WORK( TPOS ), LDT, A( I, I+KD ), LDA, ZERO, WORK( S2POS ), LDS2 );

             zhemm('Right', UPLO, PK, PN, ONE,  A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW );

             zgemm('No transpose', 'Conjugate', PK, PK, PN, ONE,  WORK( WPOS ), LDW, WORK( S2POS ), LDS2, ZERO, WORK( S1POS ), LDS1 );

             zgemm('No transpose', 'No transpose', PK, PN, PK, -HALF, WORK( S1POS ), LDS1, A( I, I+KD ), LDA, ONE,   WORK( WPOS ), LDW );


             // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
             // an update of the form:  A := A - V'*W - W'*V

             zher2k(UPLO, 'Conjugate', PN, PK, -ONE, A( I, I+KD ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA );
          } // 10

         // Copy the upper band to AB which is the band storage matrix

         for (J = N-KD+1; J <= N; J++) { // 30
            LK = min(KD, N-J) + 1;
            zcopy(LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 );
         } // 30

      } else {

          // Reduce the lower triangle of A to lower band matrix

          for (I = 1; KD < 0 ? I >= N - KD : I <= N - KD; I += KD) { // 40
             PN = N-I-KD+1;
             PK = min( N-I-KD+1, KD );

             // Compute the QR factorization of the current block

             zgeqrf(PN, KD, A( I+KD, I ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO );

             // Copy the upper portion of A into AB

             for (J = I; J <= I+PK-1; J++) { // 50
                LK = min( KD, N-J ) + 1;
                zcopy(LK, A( J, J ), 1, AB( 1, J ), 1 );
             } // 50

             zlaset('Upper', PK, PK, ZERO, ONE,  A( I+KD, I ), LDA );

             // Form the matrix T

             zlarft('Forward', 'Columnwise', PN, PK, A( I+KD, I ), LDA, TAU( I ), WORK( TPOS ), LDT );

             // Compute W:

             zgemm('No transpose', 'No transpose', PN, PK, PK, ONE, A( I+KD, I ), LDA, WORK( TPOS ), LDT, ZERO, WORK( S2POS ), LDS2 );

             zhemm('Left', UPLO, PN, PK, ONE, A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW );

             zgemm('Conjugate', 'No transpose', PK, PK, PN, ONE, WORK( S2POS ), LDS2, WORK( WPOS ), LDW, ZERO, WORK( S1POS ), LDS1 );

             zgemm('No transpose', 'No transpose', PN, PK, PK, -HALF, A( I+KD, I ), LDA, WORK( S1POS ), LDS1, ONE, WORK( WPOS ), LDW );


             // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
             // an update of the form:  A := A - V*W' - W*V'

             zher2k(UPLO, 'No transpose', PN, PK, -ONE, A( I+KD, I ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA );
             // ==================================================================
             // RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
             //  DO 45 J = I, I+PK-1
             //     LK = min( KD, N-J ) + 1
             //     CALL ZCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
// 45        CONTINUE
             // ==================================================================
          } // 40

         // Copy the lower band to AB which is the band storage matrix

         for (J = N-KD+1; J <= N; J++) { // 60
            LK = min(KD, N-J) + 1;
            zcopy(LK, A( J, J ), 1, AB( 1, J ), 1 );
         } // 60

      }

      WORK[1] = LWMIN;
      }
