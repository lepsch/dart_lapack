      SUBROUTINE DSYTRD_SY2SB( UPLO, N, KD, A, LDA, AB, LDAB, TAU,  WORK, LWORK, INFO )

      IMPLICIT NONE

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             UPLO;
      int                INFO, LDA, LDAB, LWORK, N, KD;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), AB( LDAB, * ),  TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             RONE;
      double             ZERO, ONE, HALF;
      const              RONE = 1.0D+0, ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 0.5D+0 ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY, UPPER;
      int                I, J, IINFO, LWMIN, PN, PK, LK, LDT, LDW, LDS2, LDS1, LS2, LS1, LW, LT, TPOS, WPOS, S2POS, S1POS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, DSYR2K, DSYMM, DGEMM, DCOPY, DLARFT, DGELQF, DGEQRF, DLASET
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MIN, MAX
      // ..
      // .. External Functions ..
      bool               LSAME;
      int                ILAENV2STAGE;
      // EXTERNAL LSAME, ILAENV2STAGE
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
         LWMIN = ILAENV2STAGE( 4, 'DSYTRD_SY2SB', ' ', N, KD, -1, -1 )
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
         xerbla('DSYTRD_SY2SB', -INFO );
         RETURN
      } else if ( LQUERY ) {
         WORK( 1 ) = LWMIN
         RETURN
      }

      // Quick return if possible
      // Copy the upper/lower portion of A into AB

      if ( N.LE.KD+1 ) {
          if ( UPPER ) {
              DO 100 I = 1, N
                  LK = MIN( KD+1, I )
                  dcopy(LK, A( I-LK+1, I ), 1,  AB( KD+1-LK+1, I ), 1 );
  100         CONTINUE
          } else {
              DO 110 I = 1, N
                  LK = MIN( KD+1, N-I+1 )
                  dcopy(LK, A( I, I ), 1, AB( 1, I ), 1 );
  110         CONTINUE
          ENDIF
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
      ENDIF


      // Set the workspace of the triangular matrix T to zero once such a
      // way every time T is generated the upper/lower portion will be always zero

      dlaset("A", LDT, KD, ZERO, ZERO, WORK( TPOS ), LDT );

      if ( UPPER ) {
          DO 10 I = 1, N - KD, KD
             PN = N-I-KD+1
             PK = MIN( N-I-KD+1, KD )

             // Compute the LQ factorization of the current block

             dgelqf(KD, PN, A( I, I+KD ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO );

             // Copy the upper portion of A into AB

             DO 20 J = I, I+PK-1
                LK = MIN( KD, N-J ) + 1
                dcopy(LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 );
   20        CONTINUE

             dlaset('Lower', PK, PK, ZERO, ONE,  A( I, I+KD ), LDA );

             // Form the matrix T

             dlarft('Forward', 'Rowwise', PN, PK, A( I, I+KD ), LDA, TAU( I ), WORK( TPOS ), LDT );

             // Compute W:

             dgemm('Conjugate', 'No transpose', PK, PN, PK, ONE,  WORK( TPOS ), LDT, A( I, I+KD ), LDA, ZERO, WORK( S2POS ), LDS2 );

             dsymm('Right', UPLO, PK, PN, ONE,  A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW );

             dgemm('No transpose', 'Conjugate', PK, PK, PN, ONE,  WORK( WPOS ), LDW, WORK( S2POS ), LDS2, ZERO, WORK( S1POS ), LDS1 );

             dgemm('No transpose', 'No transpose', PK, PN, PK, -HALF, WORK( S1POS ), LDS1, A( I, I+KD ), LDA, ONE,   WORK( WPOS ), LDW );


             // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
             // an update of the form:  A := A - V'*W - W'*V

             dsyr2k(UPLO, 'Conjugate', PN, PK, -ONE, A( I, I+KD ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA );
   10     CONTINUE

         // Copy the upper band to AB which is the band storage matrix

         DO 30 J = N-KD+1, N
            LK = MIN(KD, N-J) + 1
            dcopy(LK, A( J, J ), LDA, AB( KD+1, J ), LDAB-1 );
   30    CONTINUE

      } else {

          // Reduce the lower triangle of A to lower band matrix

          DO 40 I = 1, N - KD, KD
             PN = N-I-KD+1
             PK = MIN( N-I-KD+1, KD )

             // Compute the QR factorization of the current block

             dgeqrf(PN, KD, A( I+KD, I ), LDA, TAU( I ), WORK( S2POS ), LS2, IINFO );

             // Copy the upper portion of A into AB

             DO 50 J = I, I+PK-1
                LK = MIN( KD, N-J ) + 1
                dcopy(LK, A( J, J ), 1, AB( 1, J ), 1 );
   50        CONTINUE

             dlaset('Upper', PK, PK, ZERO, ONE,  A( I+KD, I ), LDA );

             // Form the matrix T

             dlarft('Forward', 'Columnwise', PN, PK, A( I+KD, I ), LDA, TAU( I ), WORK( TPOS ), LDT );

             // Compute W:

             dgemm('No transpose', 'No transpose', PN, PK, PK, ONE, A( I+KD, I ), LDA, WORK( TPOS ), LDT, ZERO, WORK( S2POS ), LDS2 );

             dsymm('Left', UPLO, PN, PK, ONE, A( I+KD, I+KD ), LDA, WORK( S2POS ), LDS2, ZERO, WORK( WPOS ), LDW );

             dgemm('Conjugate', 'No transpose', PK, PK, PN, ONE, WORK( S2POS ), LDS2, WORK( WPOS ), LDW, ZERO, WORK( S1POS ), LDS1 );

             dgemm('No transpose', 'No transpose', PN, PK, PK, -HALF, A( I+KD, I ), LDA, WORK( S1POS ), LDS1, ONE, WORK( WPOS ), LDW );


             // Update the unreduced submatrix A(i+kd:n,i+kd:n), using
             // an update of the form:  A := A - V*W' - W*V'

             dsyr2k(UPLO, 'No transpose', PN, PK, -ONE, A( I+KD, I ), LDA, WORK( WPOS ), LDW, RONE, A( I+KD, I+KD ), LDA );
             // ==================================================================
             // RESTORE A FOR COMPARISON AND CHECKING TO BE REMOVED
              // DO 45 J = I, I+PK-1
                 // LK = MIN( KD, N-J ) + 1
                 // CALL DCOPY( LK, AB( 1, J ), 1, A( J, J ), 1 )
*   45        CONTINUE
             // ==================================================================
   40     CONTINUE

         // Copy the lower band to AB which is the band storage matrix

         DO 60 J = N-KD+1, N
            LK = MIN(KD, N-J) + 1
            dcopy(LK, A( J, J ), 1, AB( 1, J ), 1 );
   60    CONTINUE

      }

      WORK( 1 ) = LWMIN
      RETURN

      // End of DSYTRD_SY2SB

      }
