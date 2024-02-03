      REAL FUNCTION SLA_GBRCOND( TRANS, N, KL, KU, AB, LDAB, AFB, LDAFB, IPIV, CMODE, C, INFO, WORK, IWORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                N, LDAB, LDAFB, INFO, KL, KU, CMODE;
      // ..
      // .. Array Arguments ..
      int                IWORK( * ), IPIV( * );
      REAL               AB( LDAB, * ), AFB( LDAFB, * ), WORK( * ), C( * )
*    ..

*  =====================================================================

      // .. Local Scalars ..
      bool               NOTRANS;
      int                KASE, I, J, KD, KE;
      REAL               AINVNM, TMP
      // ..
      // .. Local Arrays ..
      int                ISAVE( 3 );
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL SLACN2, SGBTRS, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      SLA_GBRCOND = 0.0

      INFO = 0
      NOTRANS = LSAME( TRANS, 'N' )
      if ( .NOT. NOTRANS .AND. .NOT. LSAME(TRANS, 'T') .AND. .NOT. LSAME(TRANS, 'C') ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( KL.LT.0 .OR. KL.GT.N-1 ) {
         INFO = -3
      } else if ( KU.LT.0 .OR. KU.GT.N-1 ) {
         INFO = -4
      } else if ( LDAB.LT.KL+KU+1 ) {
         INFO = -6
      } else if ( LDAFB.LT.2*KL+KU+1 ) {
         INFO = -8
      }
      if ( INFO.NE.0 ) {
         xerbla('SLA_GBRCOND', -INFO );
         RETURN
      }
      if ( N.EQ.0 ) {
         SLA_GBRCOND = 1.0
         RETURN
      }

      // Compute the equilibration matrix R such that
      // inv(R)*A*C has unit 1-norm.

      KD = KU + 1
      KE = KL + 1
      if ( NOTRANS ) {
         for (I = 1; I <= N; I++) {
            TMP = 0.0
               if ( CMODE .EQ. 1 ) {
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KD+I-J, J ) * C( J ) )
               END DO
               } else if ( CMODE .EQ. 0 ) {
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                     TMP = TMP + ABS( AB( KD+I-J, J ) )
                  END DO
               } else {
                  DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                     TMP = TMP + ABS( AB( KD+I-J, J ) / C( J ) )
                  END DO
               }
            WORK( 2*N+I ) = TMP
         END DO
      } else {
         for (I = 1; I <= N; I++) {
            TMP = 0.0
            if ( CMODE .EQ. 1 ) {
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KE-I+J, I ) * C( J ) )
               END DO
            } else if ( CMODE .EQ. 0 ) {
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KE-I+J, I ) )
               END DO
            } else {
               DO J = MAX( I-KL, 1 ), MIN( I+KU, N )
                  TMP = TMP + ABS( AB( KE-I+J, I ) / C( J ) )
               END DO
            }
            WORK( 2*N+I ) = TMP
         END DO
      }

      // Estimate the norm of inv(op(A)).

      AINVNM = 0.0

      KASE = 0
   10 CONTINUE
      slacn2(N, WORK( N+1 ), WORK, IWORK, AINVNM, KASE, ISAVE );
      if ( KASE.NE.0 ) {
         if ( KASE.EQ.2 ) {

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO

            if ( NOTRANS ) {
               sgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               sgbtrs('Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by inv(C).

            if ( CMODE .EQ. 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            } else if ( CMODE .EQ. -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }
         } else {

            // Multiply by inv(C**T).

            if ( CMODE .EQ. 1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) / C( I )
               END DO
            } else if ( CMODE .EQ. -1 ) {
               for (I = 1; I <= N; I++) {
                  WORK( I ) = WORK( I ) * C( I )
               END DO
            }

            if ( NOTRANS ) {
               sgbtrs('Transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            } else {
               sgbtrs('No transpose', N, KL, KU, 1, AFB, LDAFB, IPIV, WORK, N, INFO );
            }

            // Multiply by R.

            for (I = 1; I <= N; I++) {
               WORK( I ) = WORK( I ) * WORK( 2*N+I )
            END DO
         }
         GO TO 10
      }

      // Compute the estimate of the reciprocal condition number.

      IF( AINVNM .NE. 0.0 ) SLA_GBRCOND = ( 1.0 / AINVNM )

      RETURN

      // End of SLA_GBRCOND

      }
