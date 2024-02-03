      SUBROUTINE CHETRD_HE2HB( UPLO, N, KD, A, LDA, AB, LDAB, TAU,  WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDAB, LWORK, N, KD;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), AB( LDAB, * ),  TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               RONE
      COMPLEX            ZERO, ONE, HALF
      const              RONE = 1.0E+0, ZERO = ( 0.0E+0, 0.0E+0 ), ONE = ( 1.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                I, J, IINFO, LWMIN, PN, PK, LK, LDT, LDW, LDS2, LDS1, LS2, LS1, LW, LT, TPOS, WPOS, S2POS, S1POS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, CHER2K, CHEMM, CGEMM, CCOPY, CLARFT, CGELQF, CGEQRF, CLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      REAL               SROUNDUP_LWORK
      // EXTERNAL LSAME, ILAENV2STAGE, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Determine the minimal workspace size required
      // and test the input parameters

      INFO   = 0
      UPPER  = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      if ( N.LE.KD+1 ) {
         LWMIN = 1
      } else {
         LWMIN = ILAENV2STAGE( 4, 'CHETRD_HE2HB', '', N, KD, -1, -1 )
      }

      if ( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KD.LT.0 ) {
         INFO = -3
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -5
      } else if ( LDAB.LT.MAX( 1, KD+1 ) ) {
         INFO = -7
      } else if ( LWORK.LT.LWMIN .AND. .NOT.LQUERY ) {
         INFO = -10
      }

      if ( INFO.NE.0 ) {
         xerbla('CHETRD_HE2HB', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
         RETURN
      }

      // Quick return if possible
      // Copy the upper/lower portion of A into AB

      if ( N.LE.KD+1 ) {
          if ( UPPER ) {
              for (I = 1; I <= N; I++) { // 100
                  LK = MIN( KD+1, I )
                  ccopy(LK, A( I-LK+1, I ), 1,  AB( KD+1-LK+1, I ), 1 );
              } // 100
          } else {
              for (I = 1; I <= N; I++) { // 110
                  LK = MIN( KD+1, N-I+1 )
                  ccopy(LK, A( I, I ), 1, AB( 1, I ), 1 );
              } // 110
          }
          WORK( 1 ) = 1
          RETURN
      }

      // Determine the pointer position for the workspace

      LDT    = KD
      LDS1   = KD
      LT     = LDT*KD
      LW     = N*KD
      LS1    = LDS1*KD
      LS2    = LWMIN - LT - LW - LS1
       // LS2 = N*MAX(KD,FACTOPTNB)
      TPOS   = 1
      WPOS   = TPOS  + LT
      S1POS  = WPOS  + LW
      S2POS  = S1POS + LS1
      if ( UPPER ) {
          LDW    = KD
          LDS2   = KD
      } else {
          LDW    = N
          LDS2   = N
      }


      // Set the workspace of the triangular matrix T to zero once such a
      // way every time T is generated the upper/lower portion will be always zero

      claset("A", LDT, KD, ZERO, ZERO, WORK( TPOS ), LDT );

      if ( UPPER ) {
          DO 10 I = 1, N - KD, KD
             PN = N-I-KD+1
             PK = MIN( N-I-KD+1, KD )

             // Compute the LQ factorization of the current block

             cgelqf(KD, PN, A( I, I+KD ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO );

             // Copy the upper portion of A into AB

             DO 20 J = I, I+PK-1
                LK = MIN( KD, N-J ) + 1
                ccopy(LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 );
             } // 20

             claset('Lower', PK, PK, ZERO, ONE,  A( I, I+KD ), LDA );

             // Form the matrix T

             clarft('Forward', 'Rowwise', PN, PK, A( I, I+KD ), LDA, TAU( I ), WORK( TPOS ), LDT );

             // Compute W:

             cgemm('Conjugate', 'No transpose', PK, PN, PK, ONE,  WORK( TPOS ), LDT, A( I, I+KD ), LDA, ZERO, WORK( S2POS ), LDS2 );

             chemm('Right', UPLO, PK, PN, ONE,  A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW );

             cgemm('No transpose', 'Conjugate', PK, PK, PN, ONE,  WORK( WPOS ), LDW, WORK( S2POS ), LDS2, ZERO, WORK( S1POS ), LDS1 );

             cgemm('No transpose', 'No transpose', PK, PN, PK, -HALF, WORK( S1POS ), LDS1, A( I, I+KD ), LDA, ONE,   WORK( WPOS ), LDW );


             // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
             // an update of the form:  A := A - V'*W - W'*V

             cher2k(UPLO, 'Conjugate', PN, PK, -ONE, A( I, I+KD ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA );
          } // 10

         // Copy the upper band to AB which is the band storage matrix

         DO 30 J = N-KD+1, N
            LK = MIN(KD, N-J) + 1
            ccopy(LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 );
         } // 30

      } else {

          // Reduce the lower triangle of A to lower band matrix

          DO 40 I = 1, N - KD, KD
             PN = N-I-KD+1
             PK = MIN( N-I-KD+1, KD )

             // Compute the QR factorization of the current block

             cgeqrf(PN, KD, A( I+KD, I ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO );

             // Copy the upper portion of A into AB

             DO 50 J = I, I+PK-1
                LK = MIN( KD, N-J ) + 1
                ccopy(LK, A( J, J ), 1, AB( 1, J ), 1 );
             } // 50

             claset('Upper', PK, PK, ZERO, ONE,  A( I+KD, I ), LDA );

             // Form the matrix T

             clarft('Forward', 'Columnwise', PN, PK, A( I+KD, I ), LDA, TAU( I ), WORK( TPOS ), LDT );

             // Compute W:

             cgemm('No transpose', 'No transpose', PN, PK, PK, ONE, A( I+KD, I ), LDA, WORK( TPOS ), LDT, ZERO, WORK( S2POS ), LDS2 );

             chemm('Left', UPLO, PN, PK, ONE, A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW );

             cgemm('Conjugate', 'No transpose', PK, PK, PN, ONE, WORK( S2POS ), LDS2, WORK( WPOS ), LDW, ZERO, WORK( S1POS ), LDS1 );

             cgemm('No transpose', 'No transpose', PN, PK, PK, -HALF, A( I+KD, I ), LDA, WORK( S1POS ), LDS1, ONE, WORK( WPOS ), LDW );


             // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
             // an update of the form:  A := A - V*W' - W*V'

             cher2k(UPLO, 'No transpose', PN, PK, -ONE, A( I+KD, I ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA );
             // ==================================================================
             // RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
              // DO 45 J = I, I+PK-1
                 // LK = MIN( KD, N-J ) + 1
                 // CALL CCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
*   45        CONTINUE
             // ==================================================================
          } // 40

         // Copy the lower band to AB which is the band storage matrix

         DO 60 J = N-KD+1, N
            LK = MIN(KD, N-J) + 1
            ccopy(LK, A( J, J ), 1, AB( 1, J ), 1 );
         } // 60

      }

      WORK( 1 ) = SROUNDUP_LWORK( LWMIN )
      RETURN

      // End of CHETRD_HE2HB

      }
